def get_ref(args):
    """
    Called to get reference data from file provided in the args
    Parameters:
    args - argument values
    Functionality:
    From the reference provide, attempts to obtain reference ID and NT sequence, and AA sequence(s) if AAreport is enabled.
    For fasta formatted files, only the first entry is parsed.
    For genbank formatted files, CDS AA sequences will be pulled along with their gene ID

    Returns reference ID, nt sequence, file type and amino acid sequences encoded by the nt seq for fasta files or the CDS reported by genbank files
    """
    n = 0
    ref_id = ""
    ref_type = ""
    ref_seq = ""
    ref_orfs = {}
    if args.ref:
        ref = args.ref
        first_line = ref.readline()
        if first_line.startswith(">"):
            ref_type = "fasta"
            n += 1
            ref_id = first_line[1:].strip("\n\r").split(" ")[0]
            for line in ref:
                if line.startswith(">"):
                    n += 1
                    if n > 1:
                        break
                    ref_id = line[1:].strip("\n\r")
                elif n == 1:
                    ref_seq = ref_seq + line.strip("\n\r")
            ref_seq = ref_seq.upper()
            ref_prot = ""
            if args.AAreport == 1:
                for x in range((len(ref_seq)) // 3):
                    amino_acid = aa_call(
                        ref_seq[x * 3] + ref_seq[x * 3 + 1] + ref_seq[x * 3 + 2],
                    )
                    ref_prot = ref_prot + amino_acid
                if (len(ref_seq)) % 3 != 0:
                    ref_prot = ref_prot + "?"
                ref_orfs = [ref_id, ref_prot]
        elif first_line.upper().startswith("LOCUS"):
            ref_type = "gb"
            collect = "Null"
            orfs = {}
            reading_frames = []
            trans = 0
            nts = ""
            for line in ref:
                if collect == "Null":
                    split_line = line.strip("\n\r").split(" ")

                    if split_line[0].upper() == "VERSION":
                        ref_id = split_line[-1]
                    elif "CDS" in line:
                        collect = "CDS"
                        if "join" in line:
                            startstops = split_line[-1].strip("join()").split(",")
                            for startstop in startstops:
                                reading_frames.append(
                                    [
                                        int(startstop.split(".")[0]),
                                        int(startstop.split(".")[-1]),
                                    ],
                                )

                        else:
                            startstop = split_line[-1].split(".")
                            reading_frames.append(
                                [int(startstop[0]), int(startstop[-1])],
                            )

                    if split_line[0].upper() == "ORIGIN":
                        collect = "SEQ"
                elif collect == "CDS":
                    if "gene=" in line:
                        gene = line.split("=")[1].strip('"\n\r')
                    elif "product=" in line:
                        product = line.split("=")[1].strip('"\n\r')

                    elif "translation" in line:
                        try:
                            gene
                        except:
                            gene = "unlabeled"
                        try:
                            product
                        except:
                            product = gene

                        orf_id = product.replace(" ", "_")
                        n = 1
                        if orf_id in orfs:
                            new_id = orf_id + "." + str(n)
                            while new_id in orfs:
                                n += 1
                                new_id = orf_id + "." + str(n)
                            orf_id = new_id
                        orfs[orf_id] = {
                            "reading frames": reading_frames,
                        }
                        reading_frames = []
                        orfs[orf_id]["AAs"] = line.strip("\r\n").split('"')[1]

                        if line.strip("\r\n")[-1] != '"':
                            trans = 1
                        else:
                            orfs[orf_id]["AAs"] += "*"
                            gene = ""
                            orf_id = ""
                            collect = "Null"

                    elif trans == 1:
                        orfs[orf_id]["AAs"] = orfs[orf_id]["AAs"] + line.strip(' "\n\r')
                        if line.strip("\r\n")[-1] == '"':
                            orfs[orf_id]["AAs"] += "*"
                            trans = 0
                            gene = ""
                            orf_id = ""
                            collect = "Null"

                elif collect == "SEQ":
                    if "//" not in line:
                        ntsline = ""
                        for c in line:
                            if c.isalpha():
                                ntsline += c
                        nts += ntsline
            if args.AAreport == 1:
                for gene in orfs:
                    orfnts = ""
                    for rf in orfs[gene]["reading frames"]:
                        orfnts += nts[rf[0] - 1 : rf[1]]
                    orfs[gene]["nts"] = orfnts.upper()
            ref_orfs = orfs
            ref_seq = nts.upper()

    return (ref_id, ref_seq, ref_type, ref_orfs)


def aa_call(codon):
    """
    Called to return an amino acid enoded by the passed codon
    Parameters:
    codon - 3 nt sequence
    Functionality:
    looks up the codon in the dict and returns its value
    Returns amino acid or '?' if the codon isn't a valid codon
    """
    AADict = {
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGT": "C",
        "TGC": "C",
        "TGA": "*",
        "TGG": "W",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
        "---": "-",
    }
    AA = "?"
    if codon in AADict:
        AA = AADict[codon]

    return AA


def singlet_codon_call(nt_pos, nt, ref_nts):
    """
    Called to determine the amino acid and its position based on a single nt change relative to a referernce sequence
    Parameters:
        nt_pos - position in the reference sequence of the changed nt (1 indexed)
        nt - the new nucleotide
        ref_nts - string of the reference nts (0 indexed) encoding a peptide
    Functionality:
        calculates the position of the AA coded from the codon that the new nt is a part of and the new codon and encoded amino acid
    Returns the amino acid position (1 indexed) and the aa_call of the new codon
    """
    aa_pos = (nt_pos - 1) // 3
    aa_mod = (nt_pos - 1) % 3
    codon = ""
    try:
        if aa_mod == 0:
            codon = nt + ref_nts[nt_pos] + ref_nts[nt_pos + 1]
        elif aa_mod == 1:
            codon = ref_nts[nt_pos - 2] + nt + ref_nts[nt_pos]
        elif aa_mod == 2:
            codon = ref_nts[nt_pos - 3] + ref_nts[nt_pos - 2] + nt
    except:
        codon = "XXX"

    return (aa_pos + 1, aa_call(codon))


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
                covar_dict_dict,
                all_covar,
                file,
                args,
            )

        if file.endswith("_seqs.tsv"):
            seq_dict_dict, all_seq = collect_lines(seq_dict_dict, all_seq, file, args)

        if file.endswith("_deconv.tsv"):
            deconv_dict_dict, all_deconv = collect_lines(
                deconv_dict_dict,
                all_deconv,
                file,
                args,
            )

        if file.endswith("_pass.tsv"):
            pass_dict_dict, all_pass = collect_lines(
                pass_dict_dict,
                all_pass,
                file,
                args,
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
                            + "\t",
                        )
                    except:
                        collection_fh.write("\t\t")

                collection_fh.write("\n")
            collection_fh.close()
        else:
            print(f"No {dicts[0]} files found")
