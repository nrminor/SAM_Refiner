def collect_lines(dict_dict, all_dict, file, args):
    """
    Called to collect the lines of a file
    Parameters:
    dict_dict - dict of dict of lines from files
    all_dict - dict of all unique sequences from files
    file - name of the file to collect from
    args - argument values
    Functionality:
    Parses through the file to collect sequences and info
    Returns updated dicts
    """

    sample_line = ""
    try:
        samp = open(file)
    except:
        print("can't open " + file)
    else:
        for line in samp:
            split_line = line.strip("\n\r").split("\t")
            try:
                split_line[1]
            except:
                sample_line = split_line[0]
                dict_dict[sample_line] = {}

            else:
                if split_line[1] != "Count":
                    if float(split_line[2]) >= args.min_col_abund:
                        dict_dict[sample_line][split_line[0]] = [
                            split_line[1],
                            split_line[2],
                        ]
                        all_dict[split_line[0]] = 1
        samp.close()
        return (dict_dict, all_dict)


def collection(args):
    """
    Called to perform collection logic and printing
    Parameters:
    args - argument vlaues
    Functionality:
    Looks through directory for sample files to collect sequence info from and calls function to collect info.
    Then prints collected info to new files.
    Returns
    """

    covar_dict_dict = {}
    seq_dict_dict = {}
    deconv_dict_dict = {}
    pass_dict_dict = {}
    cr_dict_dict = {}

    all_covar = {}
    all_seq = {}
    all_deconv = {}
    all_pass = {}
    all_cr = {}

    for file in os.listdir(os.getcwd()):
        if file.endswith("_covars.tsv"):
            covar_dict_dict, all_covar = collect_lines(
                covar_dict_dict, all_covar, file, args
            )

        if file.endswith("_seqs.tsv"):
            seq_dict_dict, all_seq = collect_lines(seq_dict_dict, all_seq, file, args)

        if file.endswith("_deconv.tsv"):
            deconv_dict_dict, all_deconv = collect_lines(
                deconv_dict_dict, all_deconv, file, args
            )

        if file.endswith("_pass.tsv"):
            pass_dict_dict, all_pass = collect_lines(
                pass_dict_dict, all_pass, file, args
            )

        if file.endswith("_chim_rm.tsv"):
            cr_dict_dict, all_cr = collect_lines(cr_dict_dict, all_cr, file, args)

    dict_collection = (
        ("Covariances", covar_dict_dict, all_covar),
        ("Unique_Seqs", seq_dict_dict, all_seq),
        ("Covar_Deconv", deconv_dict_dict, all_deconv),
        ("Covar_Pass", pass_dict_dict, all_pass),
        ("Chimeras_Removed", cr_dict_dict, all_cr),
    )

    for dicts in dict_collection:
        if dicts[1]:
            if args.colID == "":
                collection_fh = open(f"Collected_{dicts[0]}.tsv", "w")
            else:
                collection_fh = open(args.colID + f"_Collected_{dicts[0]}.tsv", "w")
            sorted_seqs = sorted(dicts[2])
            collection_fh.write("\t")
            for sampline in dicts[1]:
                collection_fh.write(sampline + "\t\t")
            collection_fh.write("\nSequnces\t")
            for sampline in dicts[1]:
                collection_fh.write("Count\tAbundance\t")
            collection_fh.write("\n")
            for seq in sorted_seqs:
                collection_fh.write(seq + "\t")
                for sample in dicts[1]:
                    try:
                        collection_fh.write(
                            dicts[1][sample][seq][0]
                            + "\t"
                            + dicts[1][sample][seq][1]
                            + "\t"
                        )
                    except:
                        collection_fh.write("\t\t")

                collection_fh.write("\n")
            collection_fh.close()
        else:
            print(f"No {dicts[0]} files found")
