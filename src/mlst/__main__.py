from mlst.mlst import parse_command_line, get_input_for_mlstyper, mlstyper

def main():
    args = parse_command_line()
    input_mlstyper = get_input_for_mlstyper(args)
    mlstyper(input_mlstyper)

if __name__ == "__main__":
    main()