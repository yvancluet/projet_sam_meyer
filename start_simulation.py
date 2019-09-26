import sys
import argparse
import textwrap
import configparser
import TSC as SIM

def _main():
    parser = argparse.ArgumentParser(
      formatter_class=argparse.RawDescriptionHelpFormatter,
      epilog=textwrap.dedent('''\
         Examples:
             python3 start_simulation.py example/params.ini
             python3 start_simulation.py example/params.ini output
         '''))
    parser.add_argument("params", help="relative path to the 'params.ini' file")
    parser.add_argument("-o",
                        help="the output directory: it will be created if it doesn't exist "
                             "(default: 'first_output' or 'resume_output', "
                             "depending on whether the simulation is resumed or not)")
    args = parser.parse_args()

    try:
        ini_file = args.params
        output_dir = args.o
        SIM.start_transcribing(ini_file, output_dir)

    except configparser.NoSectionError:
        print("[Error] Please check the path to the paraneters files !")
        sys.exit(1)
    except PermissionError:
        print("[Error] Permission denied ! Please check the directory to the output files !")
        sys.exit(1)
    except (FileNotFoundError, NameError):
        print("[Error] Please check the directory to the output files !")
        sys.exit(1)


if __name__ == "__main__":
    _main()
