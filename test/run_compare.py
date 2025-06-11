import os
import subprocess


# Check if nccmp is available
def check_nccmp():
    try:
        subprocess.run(["nccmp", "--version"], capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def run_nccmp(result_file, reference_file):
    if not os.path.isfile(result_file):
        raise FileNotFoundError(f"Result file not found: {result_file}")
    if not os.path.isfile(reference_file):
        raise FileNotFoundError(f"Expected file not found: {reference_file}")

    if not check_nccmp():
        raise FileNotFoundError("nccmp not found")

    # Run nccmp with -d (data comparison) and -f (force, no user prompt)
    # Use -t for tolerance if needed (e.g., -t 1e-5 for floating-point)
    cmd = ["nccmp", "-d", "-m", "-g", "-f", "-S", result_file, reference_file]
    print(f'Testing command: {' '.join(cmd)}', flush=True)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"nccmp comparison passed: {result_file} matches {reference_file}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"nccmp comparison failed: {e.stderr}")
        return False
