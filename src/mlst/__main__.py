from mlst.mlst import parse_command_line, mlstyper

def main():
    mlstyper(parse_command_line())

if __name__ == "__main__":
    main()