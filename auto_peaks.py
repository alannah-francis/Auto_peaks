import argparse
import re

# ================================================================
# Helper: parse "Matching SO line"
# ================================================================
def parse_so_line(line):
    nums = re.findall(r'[+-]?\d+(?:\.\d+)?', line)
    if len(nums) < 3:
        return None
    try:
        so_state = int(nums[0])
        energy = float(nums[1])
    except ValueError:
        return None
    triples = []
    for i in range(2, len(nums), 3):
        if i + 2 >= len(nums):
            break
        try:
            sf_state = int(nums[i])
            spin = float(nums[i + 1])
            weight = float(nums[i + 2])
            triples.append((sf_state, spin, weight))
        except ValueError:
            continue
    return {"so_state": so_state, "energy": energy, "contributions": triples}


def format_spin(spin_value):
    return "Triplet (Spin 3)" if spin_value == 1.0 else "Singlet (Spin 1)"


# ================================================================
# Helpers for Extracted_States_Jobs.txt
# ================================================================
def parse_extracted_states(file_path):
    """
    Parse the states section from Extracted_States_Jobs.txt.
    Returns a dict keyed by SF state number -> {"JobIph": int, "Root": int}
    """
    with open(file_path, "r") as f:
        lines = f.readlines()
    state_info = {}
    states_line, jobiphs, roots = "", [], []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("State:"):
            states_line = stripped
        elif stripped.startswith("JobIph:"):
            jobiphs = [int(x) for x in re.findall(r"\d+", stripped)]
        elif stripped.startswith("Root"):
            roots = [int(x) for x in re.findall(r"\d+", stripped)]
            expected_count = len(roots)
            nums = re.findall(r"\d+", states_line)
            if len(nums) == expected_count and expected_count > 0:
                states = [int(n) for n in nums]
            else:
                # Fallback: try to split concatenated numbers in groups of 4
                merged = "".join(nums)
                if len(merged) >= 4:
                    states = [int(merged[i:i+4]) for i in range(0, len(merged), 4)]
                    states = states[:expected_count]
                else:
                    states = []
            for s, j, r in zip(states, jobiphs, roots):
                state_info[s] = {"JobIph": j, "Root": r}
    return state_info


def parse_jobiph_info(file_path):
    """
    Parse Spin/Sym info from Extracted_States_Jobs.txt.

    Priority:
      1) Look for blocks starting with:
         Specific data for JOBIPH file JOB(\d+)
         ... STATE IRREP:    (\d+)
         ... SPIN MULTIPLICITY:  (\d+)
      2) If not found for some jobs, fall back to looser patterns:
         "JobIph 4: Spin=3 Sym=2"  or "JobIph: 4 5" (presence only)
    Returns: { job_number: {"Spin": int or None, "Sym": int or None} }
    """
    jobiph_info = {}

    with open(file_path, "r") as f:
        content = f.read()

    # 1) Strong pattern: capture blocks "Specific data for JOBIPH file JOB###" ... until next block or end
    block_pattern = re.compile(
        r"Specific data for JOBIPH file JOB(\d+)(.*?)(?=Specific data for JOBIPH file JOB\d+|$)",
        re.DOTALL | re.IGNORECASE,
    )

    for m in block_pattern.finditer(content):
        jobnum = int(m.group(1))
        block = m.group(2)
        sym_m = re.search(r"STATE\s+IRREP:\s*(\d+)", block, re.IGNORECASE)
        spin_m = re.search(r"SPIN\s+MULTIPLICITY:\s*(\d+)", block, re.IGNORECASE)
        sym = int(sym_m.group(1)) if sym_m else None
        spin = int(spin_m.group(1)) if spin_m else None
        jobiph_info[jobnum] = {"Spin": spin, "Sym": sym}

    # 2) Fallback: look for explicit lines like "JobIph 4: Spin=3 Sym=2"
    #    These may appear elsewhere in the file.
    for line in content.splitlines():
        stripped = line.strip()
        m = re.match(r"JobIph\s+(\d+)\s*:\s*Spin\s*=\s*(\d+)\s+Sym\s*=\s*(\d+)", stripped, re.IGNORECASE)
        if m:
            job = int(m.group(1))
            spin = int(m.group(2))
            sym = int(m.group(3))
            jobiph_info[job] = {"Spin": spin, "Sym": sym}

    # 3) Presence-only: "JobIph: 4 5" -> ensure those jobs exist as keys (Spin/Sym unknown)
    for line in content.splitlines():
        stripped = line.strip()
        m2 = re.match(r"JobIph\s*:\s*(.*)", stripped, re.IGNORECASE)
        if m2:
            nums = [int(x) for x in re.findall(r"\d+", m2.group(1))]
            for job in nums:
                if job not in jobiph_info:
                    jobiph_info[job] = {"Spin": None, "Sym": None}

    return jobiph_info


# ================================================================
# MAIN SCRIPT
# ================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Find peaks, match to SO-RASSI states, and link SF States to Extracted_States_Jobs.txt"
    )
    parser.add_argument("peaks_file", help="Path to the intensity data file.")
    parser.add_argument("so_file", help="Path to the SO-RASSI extracted file.")
    parser.add_argument("extracted_jobs_file", help="Path to the Extracted_States_Jobs.txt file.")
    parser.add_argument("--region", type=float, nargs="*", help="Optional region center(s) for local search.")
    parser.add_argument("--window", type=int, default=10, help="Number of lines around region or global max.")
    parser.add_argument("--top_peaks", type=int, default=1, help="Number of top peaks to show.")
    parser.add_argument("--top_sf", type=int, default=1, help="Number of top spin-free states to show.")
    parser.add_argument("--output", help="Optional output file path (disables terminal printing).")
    args = parser.parse_args()

    output_lines = []
    divider = "=" * 80

    # ------------------------------------------------------------
    # LOAD PEAK DATA
    # ------------------------------------------------------------
    peaks_in_order = []
    with open(args.peaks_file, "r") as f:
        for idx, raw_line in enumerate(f, start=1):
            if not raw_line.strip():
                continue
            parts = raw_line.split()
            if len(parts) < 4:
                continue
            try:
                col1, col2, col3, col4 = parts[:4]
                intensity = float(col4)
            except ValueError:
                continue
            peaks_in_order.append((col1, col2, col3, intensity, idx))

    if not peaks_in_order:
        print("No valid peaks found.")
        return

    # ------------------------------------------------------------
    # LOAD SO-RASSI DATA
    # ------------------------------------------------------------
    so_lines = []
    with open(args.so_file, "r") as f:
        for raw in f:
            parts = raw.strip().split()
            if not parts:
                continue
            if parts[0].lstrip("-").isdigit():
                so_lines.append(raw.strip())

    def find_so_matches(final_state):
        try:
            key = int(final_state)
        except ValueError:
            return []
        matches = []
        for line in so_lines:
            parts = line.split()
            if parts and parts[0].lstrip("-").isdigit() and int(parts[0]) == key:
                matches.append(line)
        return matches

    # ------------------------------------------------------------
    # LOAD STATE & JOBIPH INFO
    # ------------------------------------------------------------
    state_info = parse_extracted_states(args.extracted_jobs_file)
    jobiph_info = parse_jobiph_info(args.extracted_jobs_file)

    # ------------------------------------------------------------
    # FUNCTION to process a list of peaks (common code for global/region)
    # ------------------------------------------------------------
    def process_peaks(peak_list, label=""):
        if label:
            output_lines.append(label)
        for col1, col2, col3, intensity, line_no in peak_list:
            output_lines.append(f"\nLine {line_no}: Transition {col1} -> {col2} : Intensity {intensity} (col3 = {col3})")
            matches = find_so_matches(col2)
            if not matches:
                output_lines.append(f"No SO match found for {col2}")
                continue

            for m in matches:
                output_lines.append(f"Matching SO line: {m}")
                parsed = parse_so_line(m)
                if parsed:
                    output_lines.append(divider)
                    output_lines.append(f"Spin-Free States of SO State {parsed['so_state']}")
                    output_lines.append(divider)
                    output_lines.append(f"Energy (au): {parsed['energy']:.6f}")

                    contributions = sorted(parsed["contributions"], key=lambda x: x[2], reverse=True)
                    for sf_state, spin, weight in contributions[: args.top_sf]:
                        output_lines.append(f"  - SF State {sf_state} | {format_spin(spin)} | Weight = {weight:.4f}")
                        if sf_state in state_info:
                            info = state_info[sf_state]
                            job_num = info["JobIph"]
                            output_lines.append(f"    JobIph: {job_num}")
                            if job_num in jobiph_info:
                                spin_mult = jobiph_info[job_num]["Spin"]
                                sym = jobiph_info[job_num]["Sym"]
                                spin_str = str(spin_mult) if spin_mult is not None else "unknown"
                                sym_str = str(sym) if sym is not None else "unknown"
                                output_lines.append(f"    Spin: {spin_str}  Sym: {sym_str}  Root: {info['Root']}")
                            else:
                                output_lines.append(f"    Root: {info['Root']}  (no JOBIPH info found)")
                        else:
                            output_lines.append("    Not found in Extracted_States_Jobs.txt")

    # ------------------------------------------------------------
    # GLOBAL MODE (no --region)
    # ------------------------------------------------------------
    if not args.region:
        global_max = max(peaks_in_order, key=lambda x: x[3])
        global_idx = global_max[4]
        global_energy = global_max[2]
        global_intensity = global_max[3]

        output_lines.append(divider)
        output_lines.append("Autopeak + State Matching Results")
        output_lines.append(divider)
        output_lines.append("")
        output_lines.append(f"Global maximum:\n  Peak #{global_idx}  Energy: {global_energy}  Intensity: {global_intensity}")
        output_lines.append("-" * 60)

        start_idx = max(0, global_idx - args.window - 1)
        end_idx = min(len(peaks_in_order), global_idx + args.window)
        local_window = peaks_in_order[start_idx:end_idx]
        local_window.sort(key=lambda x: x[3], reverse=True)
        top_peaks = local_window[: min(args.top_peaks, len(local_window))]

        process_peaks(top_peaks)

    # ------------------------------------------------------------
    # REGION MODE(s)
    # ------------------------------------------------------------
    else:
        for region_center in args.region:
            label = divider + "\n" + f"Region centered at {region_center} (±{args.window} lines)" + "\n" + divider
            # find index closest in terms of column 3 (energy)
            closest_idx = min(range(len(peaks_in_order)), key=lambda i: abs(float(peaks_in_order[i][2]) - region_center))
            start_idx = max(0, closest_idx - args.window)
            end_idx = min(len(peaks_in_order), closest_idx + args.window + 1)

            region_slice = peaks_in_order[start_idx:end_idx]
            if not region_slice:
                output_lines.append(label)
                output_lines.append(f"No peaks found near region {region_center}")
                continue

            region_slice.sort(key=lambda x: x[3], reverse=True)
            top_peaks = region_slice[: min(args.top_peaks, len(region_slice))]

            process_peaks(top_peaks, label=label)

    # ------------------------------------------------------------
    # OUTPUT RESULTS
    # ------------------------------------------------------------
    if args.output:
        with open(args.output, "w") as out:
            out.write("\n".join(output_lines))
        print(f"\nResults written to {args.output}")
    else:
        print("\n".join(output_lines))

    print("\nYIPPEE! The automated peak assignment is complete!")


if __name__ == "__main__":
    main()

