import sys
import TSC as sim

INI_file = sys.argv[1]
resume_from_path = sys.argv[2]
output_dir = sys.argv[3]

try:
	sim.resume_transcription(INI_file, resume_from_path, output_dir)
except FileNotFoundError:
	print("FILE NOT FOUND !! Please verify the path.")
except IndexError:
	print("Please specify the parameters file, file from where to resume and the output directory")

