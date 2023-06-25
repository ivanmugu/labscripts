from labscripts.mlst.mlst import parse_command_line, InputMlstyper, mlstyper

def main():
    args = parse_command_line()
    input_mlstyper = InputMlstyper(
        infile=args.infile,
        species=args.species,
        database=args.database,
        tmp_dir=args.tmp_dir,
        method_path=args.method_path,
        outdir=args.outdir,
        extented_output=args.extented_output,
        quiet=args.quiet,
        kma_matrix=args.kma_matrix,
        save_tmp=args.save_tmp,
        depth=args.depth
    )
    mlstyper(input_mlstyper)

if __name__ == "__main__":
    main()