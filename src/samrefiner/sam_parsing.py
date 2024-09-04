import itertools


def sam_line_parser(args, ref, file):
    """
    Called to pasrse through a sam file line by line and collect information for further processing
    Parameters:
    args - argument values
    ref - tuple of the information pulled from the reference file
    file - name of the sam file to be processed
    Functionality:
    opens the indicated file and runs through each line that isn't a header (starts with '@')
    and where the sequence mapped is mapped to the reference
    based on the cigar string in the line, the sequence is used to collect the nt called at each position, the indels,
    the snp variations, total number of reads mapped and the coverage
    information from each line is collected and returned
    Returns a dictionary of dictionaries for each position's nt calls, a dictionary of all the insertions,
    a list of each read id and its variant nt sequence, a dictionary of each unique variant nt sequence with total counts as values,
    the total count of reads mapped to the reference, and a dictionary for coverage of each position
    """

    samp = file[0:-4]
    nt_call_dict_dict = {}
    ins_nt_dict = {}
    reads_list = []
    col_reads = {}
    sam_read_count = 0
    sam_line_count = 0
    coverage = {}

    sam_fh = open(file)
    for line in sam_fh:
        if not line.startswith("@"):  # ignore header lines
            split_line = line.split("\t")
            if (
                ref[0].upper() == split_line[2].upper()
            ):  # check map ID matches referecne ID
                if int(split_line[4]) > 0:  # Check mapping score is positive
                    sam_line_count += 1
                    query_seq = split_line[9].upper()
                    if query_seq.strip("ATCGN-"):
                        print(
                            f"Invalid character(s) {query_seq.strip('ATCGN-')} in sequence for line {sam_line_count}, ID {split_line[0]}",
                        )
                        print("skipping")
                        continue

                    CIGAR = split_line[5]
                    cigar_num = CIGAR
                    for c in "MIDSH":
                        cigar_num = cigar_num.replace(c, "")

                    if not cigar_num.isnumeric():
                        print(
                            f"Non-standard CIGAR string {CIGAR} {cigar_num} for line {sam_line_count}, ID {split_line[0]}. Skipping",
                        )
                        continue

                    reads_count = 1
                    if args.use_count == 1:  # get the unique sequence counts
                        if "-" in split_line[0] and "=" in split_line[0]:
                            try:
                                eq_split = split_line[0].split("=")
                                dash_split = split_line[0].split("-")
                                if len(eq_split[-1]) > len(dash_split[-1]):
                                    reads_count = int(dash_split[-1])
                                else:
                                    reads_count = int(eq_split[-1])
                            except:
                                pass
                        elif "-" in split_line[0]:
                            try:
                                reads_count = int(split_line[0].split("-")[-1])
                            except:
                                pass
                        elif "=" in split_line[0]:
                            try:
                                reads_count = int(split_line[0].split("=")[-1])
                            except:
                                pass
                    sam_read_count += reads_count

                    read_start_pos = int(split_line[3])
                    readID = split_line[0]
                    run_length = 0
                    query_pos = 0
                    q_pars_pos = 0
                    mutations = []

                    for C in CIGAR:  # process sequence based on standard CIGAR line
                        if C == "M" or C == "I" or C == "D" or C == "S" or C == "H":
                            if C == "S":
                                query_pos = query_pos + run_length
                            # if C == 'H':
                            if C == "I":
                                if query_pos > 0:
                                    insert_position = q_pars_pos + read_start_pos
                                    iSeq = query_seq[query_pos : query_pos + run_length]
                                    istring = str(insert_position) + "-insert" + iSeq
                                    mutations.append(istring)
                                    if args.nt_call == 1:
                                        # add insertion to dict
                                        try:
                                            ins_nt_dict[insert_position]
                                        except:
                                            ins_nt_dict[insert_position] = {
                                                istring: reads_count,
                                            }
                                        else:
                                            try:
                                                ins_nt_dict[insert_position][
                                                    istring
                                                ] += reads_count
                                            except:
                                                ins_nt_dict[insert_position][
                                                    istring
                                                ] = reads_count

                                    query_pos = query_pos + run_length

                            elif C == "D":
                                delPos = q_pars_pos + read_start_pos
                                delnts = ref[1][delPos - 1 : delPos + run_length - 1]
                                delstring = (
                                    delnts
                                    + str(delPos)
                                    + "-"
                                    + str(delPos + run_length - 1)
                                    + "del"
                                )
                                mutations.append(delstring)
                                if args.nt_call == 1:
                                    for N in range(delPos, delPos + int(run_length)):
                                        try:
                                            nt_call_dict_dict[N]
                                        except:
                                            nt_call_dict_dict[N] = {
                                                "A": 0,
                                                "T": 0,
                                                "C": 0,
                                                "G": 0,
                                                "N": 0,
                                                "-": 0,
                                            }
                                            nt_call_dict_dict[N]["-"] = reads_count
                                        else:
                                            nt_call_dict_dict[N]["-"] += reads_count
                                q_pars_pos = q_pars_pos + run_length

                            elif C == "M":
                                offset = q_pars_pos - query_pos
                                refPos = read_start_pos + offset
                                for ntPos in range(query_pos, query_pos + run_length):
                                    if query_seq[ntPos] != ref[1][refPos + ntPos - 1]:
                                        mutations.append(
                                            ref[1][refPos + ntPos - 1]
                                            + str(refPos + ntPos)
                                            + query_seq[ntPos],
                                        )
                                    if args.nt_call == 1:
                                        if query_seq[ntPos] not in [
                                            "A",
                                            "T",
                                            "C",
                                            "G",
                                            "N",
                                            "-",
                                        ]:
                                            mut_nt = "N"
                                        else:
                                            mut_nt = query_seq[ntPos]
                                        try:
                                            nt_call_dict_dict[refPos + ntPos]
                                        except:
                                            nt_call_dict_dict[refPos + ntPos] = {
                                                "A": 0,
                                                "T": 0,
                                                "C": 0,
                                                "G": 0,
                                                "N": 0,
                                                "-": 0,
                                            }
                                            nt_call_dict_dict[refPos + ntPos][
                                                mut_nt
                                            ] = reads_count
                                        else:
                                            nt_call_dict_dict[refPos + ntPos][
                                                mut_nt
                                            ] += reads_count
                                q_pars_pos = q_pars_pos + run_length
                                query_pos = query_pos + run_length

                            run_length = 0

                        else:
                            run_length = (10 * run_length) + int(C)
                    # END CIGAR PARSE

                    seq_end_pos = read_start_pos + q_pars_pos - 1
                    if args.wgs == 1 or args.nt_call == 1:
                        for i in range(
                            read_start_pos,
                            seq_end_pos + 1,
                        ):  # update coverage
                            try:
                                coverage[i] += reads_count
                            except:
                                coverage[i] = reads_count

                    if len(mutations) == 0:  # record reference counts
                        if args.read == 1:
                            reads_list.append(
                                [readID, "Reference", read_start_pos, seq_end_pos],
                            )
                        if args.seq == 1 or args.covar == 1 or args.indel == 1:
                            try:
                                col_reads["Reference"]
                            except:
                                col_reads["Reference"] = {
                                    str(read_start_pos)
                                    + "x"
                                    + str(seq_end_pos): reads_count,
                                }
                            else:
                                try:
                                    col_reads["Reference"][
                                        str(read_start_pos) + "x" + str(seq_end_pos)
                                    ] += reads_count
                                except:
                                    col_reads["Reference"][
                                        str(read_start_pos) + "x" + str(seq_end_pos)
                                    ] = reads_count
                    else:
                        mutations = " ".join(mutations)
                        if args.read == 1:
                            reads_list.append(
                                [readID, mutations, read_start_pos, seq_end_pos],
                            )
                        if args.seq == 1 or args.covar == 1 or args.indel == 1:
                            try:
                                col_reads[mutations]
                            except:
                                col_reads[mutations] = {
                                    str(read_start_pos)
                                    + "x"
                                    + str(seq_end_pos): reads_count,
                                }
                            else:
                                try:
                                    col_reads[mutations][
                                        str(read_start_pos) + "x" + str(seq_end_pos)
                                    ] += reads_count
                                except:
                                    col_reads[mutations][
                                        str(read_start_pos) + "x" + str(seq_end_pos)
                                    ] = reads_count

    sam_fh.close()
    # END SAM LINES
    print(f"End SAM line parsing for {samp}")

    return (
        nt_call_dict_dict,
        ins_nt_dict,
        reads_list,
        col_reads,
        sam_read_count,
        coverage,
    )


def fasta_snp_call(mut, ref):
    """
    Called to determine the changed amino encoding from a single mutation event 1 indel or 1 snp based on a fasta reference
    Parameters:
    mut - string of the snp event
    ref - tuple of the information pulled from the reference file
    Functionality: for indels, determines if the change effects multiple reference codons and how. for snp, singlet_codon_call function is called to get changes
    Returns original mutation string appended with amino acid change string
    """

    if "ins" in mut:
        istring = ""
        iSeq = mut.split("insert")[1]
        insert_position = int(mut.split("-")[0])
        run_length = len(iSeq)
        if run_length % 3 == 0:
            iProt = ""
            if insert_position % 3 == 1:
                for x in range(run_length // 3):
                    AA = aa_call(iSeq[x * 3] + iSeq[x * 3 + 1] + iSeq[x * 3 + 2])
                    iProt += AA
                istring += (
                    mut + "(" + str(((insert_position - 1) // 3) + 1) + iProt + ")"
                )
            elif insert_position % 3 == 2:
                ipSeq = (
                    ref[1][insert_position - 2]
                    + iSeq
                    + ref[1][insert_position - 1 : insert_position + 1]
                )
                for x in range((run_length // 3) + 1):
                    AA = aa_call(ipSeq[x * 3] + ipSeq[x * 3 + 1] + ipSeq[x * 3 + 2])
                    iProt += AA
                istring += (
                    ref[1][insert_position - 2 : insert_position + 1]
                    + str(insert_position - 1)
                    + "-"
                    + str(insert_position + 1)
                    + ipSeq
                    + "insert("
                    + ref[3][1][(insert_position - 1) // 3]
                    + str(((insert_position - 1) // 3) + 1)
                    + iProt
                    + ")"
                )
            else:
                ipSeq = (
                    ref[1][insert_position - 3 : insert_position - 1]
                    + iSeq
                    + ref[1][insert_position - 1]
                )
                for x in range((run_length // 3) + 1):
                    AA = aa_call(ipSeq[x * 3] + ipSeq[x * 3 + 1] + ipSeq[x * 3 + 2])
                    iProt += AA
                istring += (
                    ref[1][insert_position - 3 : insert_position]
                    + str(insert_position - 2)
                    + "-"
                    + str(insert_position)
                    + ipSeq
                    + "insert("
                    + ref[3][1][(insert_position - 1) // 3]
                    + str(((insert_position - 1) // 3) + 1)
                    + iProt
                    + ")"
                )
        else:
            istring += mut + "(" + str(((insert_position - 1) // 3) + 1) + "fs)"
        return istring

    if "del" in mut:
        delPos = ""
        delnts = ""
        run_length = 0
        for c in mut.split("-")[0]:
            if c.isdigit():
                delPos += c
            else:
                run_length += 1
        delPos = int(delPos)
        delstring = ""
        delendPos = delPos + run_length - 1
        newAArefpos = (delPos - 1) // 3
        if run_length % 3 == 0:
            if (delPos) % 3 == 1:
                if run_length // 3 == 1:
                    delstring += (
                        mut
                        + "("
                        + ref[3][1][newAArefpos : newAArefpos + (run_length // 3)]
                        + str(newAArefpos + 1)
                        + "del)"
                    )
                else:
                    delstring += (
                        mut
                        + "("
                        + ref[3][1][newAArefpos : newAArefpos + (run_length // 3)]
                        + str(newAArefpos + 1)
                        + "-"
                        + str(newAArefpos + run_length // 3)
                        + "del"
                        + ")"
                    )
            else:
                if (delPos) % 3 == 2:
                    oldnts = ref[1][delPos - 2 : delendPos + 2]
                    newcodon = ref[1][delPos - 2] + ref[1][delendPos : delendPos + 2]
                    delstring += (
                        oldnts + str(delPos - 1) + "-" + str(delendPos + 2) + newcodon
                    )
                elif (delPos) % 3 == 0:
                    oldnts = ref[1][delPos - 3 : delendPos + 1]
                    newcodon = ref[1][delPos - 3 : delPos - 1] + ref[1][delendPos]
                    delstring += (
                        oldnts + str(delPos - 2) + "-" + str(delendPos + 1) + newcodon
                    )
                delstring += (
                    "del("
                    + ref[3][1][newAArefpos : newAArefpos + (run_length // 3) + 1]
                    + str(newAArefpos + 1)
                    + "-"
                    + str(newAArefpos + 1 + run_length // 3)
                    + aa_call(newcodon)
                    + "del)"
                )
        else:
            delstring += mut + "(" + str(newAArefpos + 1) + "fs)"
        return delstring
    AAinfo = singlet_codon_call(int(mut[1:-1]), mut[-1], ref[1])
    AAstring = "(" + ref[3][1][AAinfo[0] - 1] + str(AAinfo[0]) + AAinfo[1] + ")"
    return mut + AAstring


def gb_snp_call(mut, ref):
    """
    Called to determine the changed amino encoding from a single mutation event 1 indel or 1 snp based on a gb reference
    Parameters:
    mut - string of the snp event
    ref - tuple of the information pulled from the reference file
    Functionality: Determines amino acid changes for each orf effected.  For indels, determines if the change effects multiple reference codons and how.
    for snp, singlet_codon_call function is called to get changes
    Returns original mutation string appended with amino acid change string
    """

    if "ins" in mut:
        iSeq = mut.split("insert")[1]
        insert_position = int(mut.split("-")[0])
        run_length = len(iSeq)
        iProt = ""
        for orf in ref[3]:
            orflength = 0
            for rf in ref[3][orf]["reading frames"]:
                if insert_position >= rf[0] and insert_position <= rf[1]:
                    iProt += "|(" + orf + ":"
                    orfPos = 1 + insert_position - rf[0] + orflength
                    AA = ""
                    if run_length % 3 == 0:
                        if orfPos % 3 == 1:
                            for x in range(run_length // 3):
                                AA += aa_call(
                                    iSeq[x * 3] + iSeq[x * 3 + 1] + iSeq[x * 3 + 2],
                                )
                            iProt += (
                                str(orfPos)
                                + "insert"
                                + iSeq
                                + "("
                                + str(((orfPos - 1) // 3) + 1)
                                + AA
                            )
                        else:
                            if orfPos % 3 == 2:
                                ipSeq = (
                                    ref[1][insert_position - 2]
                                    + iSeq
                                    + ref[1][insert_position - 1 : insert_position + 1]
                                )
                                for x in range((run_length // 3) + 1):
                                    AA += aa_call(
                                        ipSeq[x * 3]
                                        + ipSeq[x * 3 + 1]
                                        + ipSeq[x * 3 + 2],
                                    )
                                iProt += (
                                    ref[1][insert_position - 2 : insert_position + 1]
                                    + str(orfPos - 1)
                                    + "-"
                                    + str(orfPos + 1)
                                    + ipSeq
                                )
                            else:
                                ipSeq = (
                                    ref[1][insert_position - 3 : insert_position - 1]
                                    + iSeq
                                    + ref[1][insert_position - 1]
                                )
                                for x in range((run_length // 3) + 1):
                                    AA += aa_call(
                                        ipSeq[x * 3]
                                        + ipSeq[x * 3 + 1]
                                        + ipSeq[x * 3 + 2],
                                    )
                                iProt += (
                                    ref[1][insert_position - 3 : insert_position]
                                    + str(orfPos - 2)
                                    + "-"
                                    + str(orfPos)
                                    + ipSeq
                                )
                            iProt += (
                                "insert("
                                + ref[3][orf]["AAs"][(orfPos - 1) // 3]
                                + str(((orfPos - 1) // 3) + 1)
                                + AA
                            )
                    else:
                        iProt += (
                            str(orfPos)
                            + "insert"
                            + iSeq
                            + "("
                            + str(((orfPos - 1) // 3) + 1)
                            + "fs"
                        )
                    iProt += "))"
                orflength += rf[1] - rf[0] + 1
        return mut + iProt

    if "del" in mut:
        delPos = ""
        delnts = ""
        run_length = 0
        for c in mut.split("-")[0]:
            if c.isdigit():
                delPos += c
            else:
                run_length += 1
        delPos = int(delPos)
        delProt = ""
        delendPos = delPos + run_length - 1
        for orf in ref[3]:
            orflength = 0
            for rf in ref[3][orf]["reading frames"]:
                if delPos > rf[1] or delendPos < rf[0]:
                    pass
                elif delPos >= rf[0] and delendPos <= rf[1]:
                    delProt += "|(" + orf + ":"
                    orfPos = 1 + delPos - rf[0] + orflength
                    orfendPos = orfPos + run_length - 1
                    startcodon = ((orfPos - 1) // 3) + 1
                    if run_length % 3 == 0:
                        endcodon = ((orfendPos - 1) // 3) + 1
                        if (orfPos) % 3 == 1:
                            delProt += (
                                delnts
                                + str(orfPos)
                                + "-"
                                + str(orfendPos)
                                + "del("
                                + ref[3][orf]["AAs"][startcodon - 1 : endcodon]
                                + str(startcodon)
                            )
                            if run_length // 3 == 1:
                                delProt += "del"
                            else:
                                delProt += "-" + str(endcodon) + "del"
                        else:
                            delAA = ""
                            if (orfPos) % 3 == 2:
                                oldnts = ref[1][delPos - 2 : delendPos + 2]
                                newcodon = (
                                    ref[1][delPos - 2]
                                    + ref[1][delendPos : delendPos + 2]
                                )
                                delProt += (
                                    oldnts
                                    + str(orfPos - 1)
                                    + "-"
                                    + str(orfendPos + 2)
                                    + newcodon
                                )
                            else:
                                oldnts = ref[1][delPos - 3 : delendPos + 1]
                                newcodon = (
                                    ref[1][delPos - 3 : delPos - 1] + ref[1][delendPos]
                                )
                                delProt += (
                                    oldnts
                                    + str(orfPos - 2)
                                    + "-"
                                    + str(orfendPos + 1)
                                    + newcodon
                                )
                            delProt += (
                                "del("
                                + ref[3][orf]["AAs"][startcodon - 1 : endcodon]
                                + str(startcodon)
                                + "-"
                                + str(endcodon)
                                + aa_call(newcodon)
                                + "del"
                            )
                    else:
                        delProt += (
                            delnts
                            + str(orfPos)
                            + "-"
                            + str(orfendPos)
                            + "del("
                            + str(startcodon)
                            + "fs"
                        )
                    delProt += "))"
                elif delendPos >= rf[0] and delPos <= rf[0]:
                    delProt += "|(" + orf
                    if orflength == 0:
                        delProt += ":Start_disrupted)"
                    else:
                        delProt += ":fs/Splicing_disrupted)"
                else:
                    delProt += "|(" + orf + ":fs/Splicing/Termination_disrupted)"
                orflength += rf[1] - rf[0] + 1
        return mut + delProt
    PM = ""
    PMPos = int(mut[1:-1])
    for orf in ref[3]:
        orflength = 0
        for rf in ref[3][orf]["reading frames"]:
            if PMPos >= rf[0] and PMPos <= rf[1]:
                ORFPos = PMPos - rf[0] + 1 + orflength
                AAinfo = singlet_codon_call(ORFPos, mut[-1], (ref[3][orf]["nts"]))
                PM += (
                    "|("
                    + orf
                    + ":"
                    + mut[0]
                    + str(ORFPos)
                    + mut[-1]
                    + "("
                    + ref[3][orf]["AAs"][AAinfo[0] - 1]
                    + str(AAinfo[0])
                    + AAinfo[1]
                    + "))"
                )
            orflength += rf[1] - rf[0] + 1
    return mut + PM


def get_combos(qlist, clen):
    """
    Called to get combinations of cosegregating polymorphisms
    Parameters:
    qlist - list of polymorphisms
    clen - number of polymorphisms to report together

    Functionality:
        uses itertools.combinations to get sets of combinations for each length
    Returns list of combinations
    """
    combos = []
    if clen == 0 or clen > len(qlist):
        clen = len(qlist)
    for N in range(1, clen + 1):
        for comb in itertools.combinations(qlist, N):
            combos.append(" ".join(comb))
    return combos


def print_indels(samp, sam_read_count, indel_dict, coverage, args):
    """
    Called to print the indel output
    Parameters:
    samp - name of the sam being processed
    sam_read_count - number of total reads mapped in the sam file
    indel_reads - dict of indels
    coverage - dict of coverage at reference positions
    args - arguement values

    Functionality:
        sorts indels and gets abundance for wgs mode or not, if any indels pass abundance check, opens file and writes lines
    Returns nothing
    """
    sorted_indels = sorted(indel_dict, key=indel_dict.__getitem__, reverse=True)
    indels_to_write = []
    for key in sorted_indels:
        if indel_dict[key] >= args.min_count:
            if args.wgs == 0:
                abund = indel_dict[key] / sam_read_count
            elif args.wgs == 1:
                indelPos = ""
                for c in key.strip("ATCGN"):
                    if c.isdigit():
                        indelPos += c
                    else:
                        break
                abund = indel_dict[key] / coverage[int(indelPos)]
            if abund >= args.min_samp_abund:
                indels_to_write.append(f"{key}\t{indel_dict[key]}\t{abund:.3f}\n")
        else:
            break
    if len(indels_to_write) > 0:
        indel_fh = open(samp + "_indels.tsv", "w")
        indel_fh.write(samp + "(" + str(sam_read_count) + ")\n")
        indel_fh.write("Indel\tCount\tAbundance\n")
        for indel_entry in indels_to_write:
            indel_fh.write(indel_entry)
        indel_fh.close()


def print_covars(samp, sam_read_count, col_reads, coverage, args, aa_centered):
    """
    Called to print the covars output
    Parameters:
    samp - name of the sam being processed
    sam_read_count - number of total reads mapped in the sam file
    col_reads - dict of collapsed sequences with amino acid seq if applicable
    coverage - dict of coverage at reference positions
    args - arguement values
    aa_centered- dict of coverage at posistion of polymoprhisms in the aa centered form

    Functionality:
        parses each unique variant sequence, normal and aa centered, if applicable, to get the cosegregating variations and tiled coverage if enabled,
        if there are combos that meet the count threshold output files are opened and lines written
    Returns nothing
    """
    ntcombinations = {}
    AAcombinations = {}
    tiles = {}
    for sequence in col_reads:
        ntsingles = []
        AAsingles = []

        if args.wgs == 1 and sequence == "Reference":
            continue
        try:
            ntsingles = col_reads[sequence]["AA_sequence"].split(" ")
        except:
            ntsingles = sequence.split(" ")
        if aa_centered:
            AAsingles = col_reads[sequence]["AA_centered"].split(" ")
        sequence_count = 0
        for value in col_reads[sequence].values():
            if isinstance(value, int):
                sequence_count += value

        if args.wgs == 1 and args.covar_tile_coverage == 1:
            for start_end in col_reads[SNP_sequence]:
                if "x" in start_end:
                    start_pos, end_pos = start_end.split("x")
                    start_pos = int(start_pos)
                    end_pos = int(end_pos)
                    try:
                        tiles[start_pos]
                    except:
                        tiles[start_pos] = {end_pos: col_reads[SNP_sequence][start_end]}
                    else:
                        try:
                            tiles[start_pos][end_pos] += col_reads[SNP_sequence][
                                start_end
                            ]
                        except:
                            tiles[start_pos][end_pos] = col_reads[SNP_sequence][
                                start_end
                            ]
        if ntsingles and len(ntsingles) <= args.max_dist:
            for combo in get_combos(ntsingles, args.max_covar):
                if combo not in ntcombinations:
                    ntcombinations[combo] = sequence_count
                else:
                    ntcombinations[combo] += sequence_count

        if AAsingles and len(AAsingles) <= args.max_dist:
            for combo in get_combos(AAsingles, args.max_covar):
                if combo not in AAcombinations:
                    AAcombinations[combo] = sequence_count
                else:
                    AAcombinations[combo] += sequence_count
    for combinations, samp_tag in ((ntcombinations, ""), (AAcombinations, "_AA")):
        try:
            max_count = max(combinations.values())
        except:
            pass
        else:
            if max_count >= args.min_count:
                covar_fh = open(samp + samp_tag + "_covars.tsv", "w")
                covar_fh.write(samp + "(" + str(sam_read_count) + ")\n")
                covar_fh.write("Co-Variants\tCount\tAbundance\n")
                sortedcombos = sorted(
                    combinations,
                    key=combinations.__getitem__,
                    reverse=True,
                )
                for key in sortedcombos:
                    if combinations[key] >= args.min_count:
                        if args.wgs == 0:
                            abund = combinations[key] / sam_read_count
                        elif args.wgs == 1:
                            splitcombos = key.split()
                            if not samp_tag:
                                startcovPos = ""
                                for c in splitcombos[0].strip("ATGC"):
                                    if c.isdigit():
                                        startcovPos += c
                                    else:
                                        break
                                endcovPos = ""
                                split_mut = splitcombos[-1].split("(")[0]
                                if "del" in split_mut:
                                    pos_string = split_mut.strip("ATGCN").split("-")[1]
                                else:
                                    pos_string = split_mut.strip("ATGCN")
                                for c in pos_string:
                                    if c.isdigit():
                                        endcovPos += c
                                    else:
                                        break
                                if not endcovPos:
                                    print(key)
                                    print(split_mut)
                                    print(pos_string)
                            else:
                                startcovPos = aa_centered[splitcombos[0]][0]
                                endcovPos = aa_centered[splitcombos[-1]][-1]
                            if startcovPos == endcovPos:
                                abund = combinations[key] / coverage[int(startcovPos)]

                            elif args.covar_tile_coverage == 0:
                                coveragevals = []
                                for i in range(int(startcovPos), int(endcovPos) + 1):
                                    coveragevals.append(coverage[i])
                                    abund = combinations[key] / min(coveragevals)
                            elif args.covar_tile_coverage == 1:
                                coverageval = 0
                                for tile_start in tiles:
                                    if int(startcovPos) >= int(tile_start):
                                        for tile_end in tiles[tile_start]:
                                            if int(endcovPos) <= int(tile_end):
                                                coverageval += tiles[tile_start][
                                                    tile_end
                                                ]
                                abund = combinations[key] / coverageval
                        if abund >= args.min_samp_abund:
                            covar_fh.write(f"{key}\t{combinations[key]}\t{abund:.3f}\n")
                    else:
                        break

                covar_fh.close()


def print_covars(samp, sam_read_count, col_reads, coverage, args, aa_centered):
    """
    Called to print the covars output
    Parameters:
    samp - name of the sam being processed
    sam_read_count - number of total reads mapped in the sam file
    col_reads - dict of collapsed sequences with amino acid seq if applicable
    coverage - dict of coverage at reference positions
    args - arguement values
    aa_centered- dict of coverage at posistion of polymoprhisms in the aa centered form

    Functionality:
        parses each unique variant sequence, normal and aa centered, if applicable, to get the cosegregating variations and tiled coverage if enabled,
        if there are combos that meet the count threshold output files are opened and lines written
    Returns nothing
    """
    ntcombinations = {}
    AAcombinations = {}
    tiles = {}
    for sequence in col_reads:
        ntsingles = []
        AAsingles = []

        if args.wgs == 1 and sequence == "Reference":
            continue
        try:
            ntsingles = col_reads[sequence]["AA_sequence"].split(" ")
        except:
            ntsingles = sequence.split(" ")
        if aa_centered:
            AAsingles = col_reads[sequence]["AA_centered"].split(" ")
        sequence_count = 0
        for value in col_reads[sequence].values():
            if isinstance(value, int):
                sequence_count += value

        if args.wgs == 1 and args.covar_tile_coverage == 1:
            for start_end in col_reads[SNP_sequence]:
                if "x" in start_end:
                    start_pos, end_pos = start_end.split("x")
                    start_pos = int(start_pos)
                    end_pos = int(end_pos)
                    try:
                        tiles[start_pos]
                    except:
                        tiles[start_pos] = {end_pos: col_reads[SNP_sequence][start_end]}
                    else:
                        try:
                            tiles[start_pos][end_pos] += col_reads[SNP_sequence][
                                start_end
                            ]
                        except:
                            tiles[start_pos][end_pos] = col_reads[SNP_sequence][
                                start_end
                            ]
        if ntsingles and len(ntsingles) <= args.max_dist:
            for combo in get_combos(ntsingles, args.max_covar):
                if combo not in ntcombinations:
                    ntcombinations[combo] = sequence_count
                else:
                    ntcombinations[combo] += sequence_count

        if AAsingles and len(AAsingles) <= args.max_dist:
            for combo in get_combos(AAsingles, args.max_covar):
                if combo not in AAcombinations:
                    AAcombinations[combo] = sequence_count
                else:
                    AAcombinations[combo] += sequence_count
    for combinations, samp_tag in ((ntcombinations, ""), (AAcombinations, "_AA")):
        try:
            max_count = max(combinations.values())
        except:
            pass
        else:
            if max_count >= args.min_count:
                covar_fh = open(samp + samp_tag + "_covars.tsv", "w")
                covar_fh.write(samp + "(" + str(sam_read_count) + ")\n")
                covar_fh.write("Co-Variants\tCount\tAbundance\n")
                sortedcombos = sorted(
                    combinations,
                    key=combinations.__getitem__,
                    reverse=True,
                )
                for key in sortedcombos:
                    if combinations[key] >= args.min_count:
                        if args.wgs == 0:
                            abund = combinations[key] / sam_read_count
                        elif args.wgs == 1:
                            splitcombos = key.split()
                            if not samp_tag:
                                startcovPos = ""
                                for c in splitcombos[0].strip("ATGC"):
                                    if c.isdigit():
                                        startcovPos += c
                                    else:
                                        break
                                endcovPos = ""
                                split_mut = splitcombos[-1].split("(")[0]
                                if "del" in split_mut:
                                    pos_string = split_mut.strip("ATGCN").split("-")[1]
                                else:
                                    pos_string = split_mut.strip("ATGCN")
                                for c in pos_string:
                                    if c.isdigit():
                                        endcovPos += c
                                    else:
                                        break
                                if not endcovPos:
                                    print(key)
                                    print(split_mut)
                                    print(pos_string)
                            else:
                                startcovPos = aa_centered[splitcombos[0]][0]
                                endcovPos = aa_centered[splitcombos[-1]][-1]
                            if startcovPos == endcovPos:
                                abund = combinations[key] / coverage[int(startcovPos)]

                            elif args.covar_tile_coverage == 0:
                                coveragevals = []
                                for i in range(int(startcovPos), int(endcovPos) + 1):
                                    coveragevals.append(coverage[i])
                                    abund = combinations[key] / min(coveragevals)
                            elif args.covar_tile_coverage == 1:
                                coverageval = 0
                                for tile_start in tiles:
                                    if int(startcovPos) >= int(tile_start):
                                        for tile_end in tiles[tile_start]:
                                            if int(endcovPos) <= int(tile_end):
                                                coverageval += tiles[tile_start][
                                                    tile_end
                                                ]
                                abund = combinations[key] / coverageval
                        if abund >= args.min_samp_abund:
                            covar_fh.write(f"{key}\t{combinations[key]}\t{abund:.3f}\n")
                    else:
                        break

                covar_fh.close()


def get_nt_indels(col_reads):
    """
    Called to collect indels if amino acid reporting is disabled
    Parameters:
    col_reads - dict of unique variant sequences
    Functionality: goes through unique sequences looking for indels and adds them to a new dict
    Returns the new indel dict
    """
    indel_dict = {}
    for SNP_sequence in col_reads:
        if "del" in SNP_sequence or "insert" in SNP_sequence:
            for mut in SNP_sequence.split(" "):
                if "del" in mut or "ins" in mut:
                    try:
                        indel_dict[mut] += sum(col_reads[SNP_sequence].values())
                    except:
                        indel_dict[mut] = sum(col_reads[SNP_sequence].values())
    return indel_dict


def fa_sam_parse(args, ref, file):
    """
    Called to handle main sam info parsing logic when using a fasta reference
    Parameters:
    args - argument values
    ref - tuple of referecne information
    file - name of sam file
    Functionality:
    calls sam_line_parser to get info from sam, then processes the information to get amino acid change information if applicable.
    prints outputs, mainly by function calls.  nt call output is the exception
    Returns nothing
    """

    samp = file[0:-4]
    print(f"Starting {samp} processing")
    nt_call_dict_dict, ins_nt_dict, reads_list, col_reads, sam_read_count, coverage = (
        sam_line_parser(args, ref, file)
    )

    if sam_read_count == 0:
        print(f"No Reads for {samp}")
    else:
        indel_dict = {}
        nt_to_AA_dict = {}

        if args.AAreport == 1 and (
            args.seq == 1 or args.read == 1 or args.covar == 1 or args.indel == 1
        ):
            for SNP_sequence in col_reads:
                if SNP_sequence == "Reference":
                    pass
                elif args.AAcodonasMNP == 1:
                    try:
                        col_reads[SNP_sequence]["AA_sequence"]
                    except:
                        MNPs = []
                        curMNP = ""
                        last_codon = -1
                        fshift = 0
                        for mut in SNP_sequence.split(
                            " ",
                        ):  # collect together mutations that affect the same codon(s)
                            startPos = int(mut.split("-")[0].strip("ATCGN"))
                            endPos = startPos
                            if "del" in mut:
                                endPos = int(mut.split("-")[1].strip("del"))

                            startcodon = (((startPos) - 1) // 3) + 1
                            endcodon = ((endPos - 1) // 3) + 1
                            if curMNP:
                                if "fs" in mut:
                                    mutshift = 0
                                    if "insert" in mut:
                                        mutshift += len(
                                            mut.split("(")[0].split("insert")[1],
                                        )
                                    else:
                                        mutshift -= len(
                                            mut.split("-")[0].strip("0123456789"),
                                        )

                                    if (fshift % 3) == 0:
                                        if startcodon == last_codon:
                                            curMNP.append([mut, startPos])
                                            fshift += mutshift
                                        else:
                                            MNPs.append(curMNP)
                                            curMNP = [[mut, startPos]]
                                            fshift = mutshift

                                    elif startcodon < last_codon + 5:
                                        curMNP.append([mut, startPos])
                                        fshift += mutshift
                                    else:
                                        ## TODO splitfsMNP() make process
                                        fsfreeMNP = []
                                        for SNP in curMNP:
                                            if "fs" in SNP:
                                                if fsfreeMNP:
                                                    MNPs.append(fsfreeMNP)
                                                    fsfreeMNP = []
                                                MNPs.append([SNP])
                                            else:
                                                fsfreeMNP.append(SNP)
                                        if fsfreeMNP:
                                            MNPs.append(fsfreeMNP)
                                        curMNP = [[mut, startPos]]
                                        fshift = mutshift
                                elif startcodon == last_codon:
                                    curMNP.append([mut, startPos])
                                else:
                                    if (fshift % 3) == 0:
                                        MNPs.append(curMNP)
                                    else:
                                        ## splitfsMNP() make process
                                        fsfreeMNP = []
                                        for SNP in curMNP:
                                            if "fs" in SNP:
                                                if fsfreeMNP:
                                                    MNPs.append(fsfreeMNP)
                                                    fsfreeMNP = []
                                                MNPs.append([SNP])
                                            else:
                                                fsfreeMNP.append(SNP)
                                        if fsfreeMNP:
                                            MNPs.append(fsfreeMNP)
                                    curMNP = [[mut, startPos]]
                                    fshift = 0
                            else:
                                curMNP = [[mut, startPos]]
                            last_codon = endcodon
                        if curMNP:
                            MNPs.append(curMNP)

                        if MNPs:
                            combinedmut = []
                            for entry in MNPs:
                                if (
                                    len(entry) > 1
                                ):  # re-calcs sequence for codons affected by multiple changes
                                    nt_to_AA_key = []
                                    for snpmut in entry:
                                        nt_to_AA_key.append(snpmut[0])
                                    nt_to_AA_key = " ".join(nt_to_AA_key)
                                    try:
                                        newmut = nt_to_AA_dict[nt_to_AA_key]
                                    except:
                                        MNPreport = ""
                                        mutntseq = ""
                                        orfstartpos = entry[0][1]
                                        orfendpos = entry[-1][1]
                                        if "del" in entry[-1][0]:
                                            orfendpos += (
                                                len(
                                                    entry[-1][0]
                                                    .split("-")[0]
                                                    .strip("0123456789"),
                                                )
                                                - 1
                                            )

                                        for i in range(len(entry)):
                                            curPM = entry[i][0]
                                            if i > 0:
                                                if entry[i][1] > entry[i - 1][1]:
                                                    if "insert" in entry[i - 1][0]:
                                                        mutntseq += ref[1][
                                                            entry[i - 1][1] - 1 : entry[
                                                                i
                                                            ][1]
                                                            - 1
                                                        ]
                                                    elif "del" in entry[i - 1][0]:
                                                        mutntseq += ref[1][
                                                            entry[i - 1][1]
                                                            + len(
                                                                entry[i - 1][0]
                                                                .split("-")[0]
                                                                .strip("0123456789"),
                                                            )
                                                            - 1 : entry[i][1] - 1
                                                        ]
                                                    else:
                                                        mutntseq += ref[1][
                                                            entry[i - 1][1] : entry[i][
                                                                1
                                                            ]
                                                            - 1
                                                        ]
                                            if "insert" in curPM:
                                                mutntseq += curPM.split("(")[0].split(
                                                    "insert",
                                                )[1]
                                            elif "del" in curPM:
                                                pass
                                            else:
                                                mutntseq += curPM[-1]
                                        startcodonpos = ((orfstartpos - 1) // 3) + 1
                                        endcodonpos = ((orfendpos - 1) // 3) + 1
                                        startmod = orfstartpos % 3
                                        endmod = orfendpos % 3
                                        wtorfstartpos = orfstartpos
                                        wtorfendpos = orfendpos
                                        if startmod == 2:
                                            wtorfstartpos = orfstartpos - 1
                                            mutntseq = (
                                                ref[1][orfstartpos - 2] + mutntseq
                                            )
                                        elif startmod == 0:
                                            wtorfstartpos = orfstartpos - 2
                                            mutntseq = (
                                                ref[1][
                                                    orfstartpos - 3 : orfstartpos - 1
                                                ]
                                                + mutntseq
                                            )
                                        if endmod == 2:
                                            wtorfendpos = orfendpos + 1
                                            mutntseq += ref[1][orfendpos]
                                        elif endmod == 1:
                                            wtorfendpos = orfendpos + 2
                                            mutntseq += ref[1][
                                                orfendpos : orfendpos + 2
                                            ]
                                        wtntseq = ref[1][
                                            wtorfstartpos - 1 : wtorfendpos
                                        ]
                                        wtAAseq = ref[3][1][
                                            startcodonpos - 1 : endcodonpos
                                        ]
                                        mutAAseq = ""
                                        for i in range(len(mutntseq) // 3):
                                            mutAAseq += aa_call(
                                                mutntseq[i * 3 : (i * 3) + 3],
                                            )
                                        if len(mutntseq) % 3 != 0:
                                            mutAAseq += "fs"
                                        if len(mutAAseq) < len(wtAAseq):
                                            mutAAseq += "del"
                                            mutntseq += "del"
                                        elif len(mutAAseq) > len(wtAAseq):
                                            mutAAseq += "insert"
                                            mutntseq += "insert"
                                        if startcodonpos == endcodonpos:
                                            newmutstring = f"{wtntseq}{wtorfstartpos}-{wtorfendpos}{mutntseq}({wtAAseq}{startcodonpos}{mutAAseq})"
                                        else:
                                            newmutstring = f"{wtntseq}{wtorfstartpos}-{wtorfendpos}{mutntseq}({wtAAseq}{startcodonpos}-{endcodonpos}{mutAAseq})"

                                        nt_to_AA_dict[nt_to_AA_key] = newmutstring
                                        newmut = newmutstring
                                elif entry[0][0] in nt_to_AA_dict:
                                    newmut = nt_to_AA_dict[entry[0][0]]
                                else:
                                    mutstring = fasta_snp_call(entry[0][0], ref)
                                    nt_to_AA_dict[entry[0][0]] = mutstring
                                    newmut = mutstring
                                combinedmut.append(newmut)
                                if args.indel == 1:
                                    if "del" in newmut or "insert" in newmut:
                                        try:
                                            indel_dict[newmut] += sum(
                                                col_reads[SNP_sequence].values(),
                                            )
                                        except:
                                            indel_dict[newmut] = sum(
                                                col_reads[SNP_sequence].values(),
                                            )
                        col_reads[SNP_sequence]["AA_sequence"] = " ".join(combinedmut)
                else:
                    try:
                        col_reads[SNP_sequence]["AA_sequence"]
                    except:
                        mutations = []
                        for mut in SNP_sequence.split(" "):
                            try:
                                newmut = nt_to_AA_dict[mut]
                            except:
                                mutstring = fasta_snp_call(mut, ref)
                                nt_to_AA_dict[mut] = mutstring
                                newmut = mutstring
                            mutations.append(newmut)
                            if args.indel == 1:
                                if "del" in mut or "insert" in mut:
                                    try:
                                        indel_dict[newmut] += sum(
                                            col_reads[SNP_sequence].values(),
                                        )
                                    except:
                                        indel_dict[newmut] = sum(
                                            col_reads[SNP_sequence].values(),
                                        )
                        col_reads[SNP_sequence]["AA_sequence"] = " ".join(mutations)

        elif args.indel == 1:
            indel_dict = get_nt_indels(col_reads)

        if args.read == 1:
            print_reads(samp, reads_list, col_reads, args)

        if args.seq == 1:
            print_unique_seq(samp, sam_read_count, col_reads, coverage, args)

        if args.indel == 1 and len(indel_dict) > 0:
            print_indels(samp, sam_read_count, indel_dict, coverage, args)

        if args.covar == 1:
            print_covars(samp, sam_read_count, col_reads, coverage, args, {})

        if args.nt_call == 1:
            ntcall_lines = {
                "line": {},
                "variant": {},
            }
            ntcall_fh = open(samp + "_nt_calls.tsv", "w")
            ntcall_fh.write(samp + "(" + str(sam_read_count) + ")\n")
            if args.ntvar == 1:
                ntcallv_fh = open(samp + "_nt_calls_varonly.tsv", "w")
                ntcallv_fh.write(samp + "(" + str(sam_read_count) + ")\n")
            sorted_Pos = sorted(nt_call_dict_dict)
            if args.AAreport == 1:
                ntcall_fh.write(
                    "Position\tref NT\tAA Pos\tref AA\tA\tT\tC\tG\t-\tN\tTotal\tPrimary NT\tCounts\tAbundance\tPrimary Seq AA\tsingle nt AA\tSecondary NT\tCounts\tAbundance\tAA\tTertiary NT\tCounts\tAbundance\tAA\n",
                )
                if args.ntvar == 1:
                    ntcallv_fh.write(
                        "Position\tref NT\tAA Pos\tref AA\tA\tT\tC\tG\t-\tN\tTotal\tPrimary NT\tCounts\tAbundance\tPrimary Seq AA\tsingle nt AA\tSecondary NT\tCounts\tAbundance\tAA\tTertiary NT\tCounts\tAbundance\tAA\n",
                    )
                consensus = {}
                for Pos in sorted_Pos:
                    try:
                        total = coverage[Pos]
                    except:
                        total = 0
                        print(f"coverage of position {Pos} not found")

                    if (
                        total >= (sam_read_count * args.ntabund) or args.wgs == 1
                    ) and total >= args.ntcover:
                        try:
                            ins_nt_dict[Pos]
                        except:
                            pass
                        else:
                            i = 1
                            for insertion in ins_nt_dict[Pos]:
                                insert_position = Pos + (i / 1000)
                                i_nts = insertion.split("insert")[1]
                                try:
                                    i_AAs = (
                                        nt_to_AA_dict[insertion]
                                        .split("(")[1]
                                        .strip(")")
                                    )
                                except:
                                    i_AAs = (
                                        fasta_snp_call(insertion, ref)
                                        .split("(")[1]
                                        .strip(")")
                                    )
                                AA_pos = ((Pos - 1) // 3) + 1
                                iabund = ins_nt_dict[Pos][insertion] / total
                                ntcall_lines["line"][insert_position] = (
                                    f"{Pos}\t-\t{AA_pos}\t-\t\t\t\t\t\t\t{total}\t{i_nts}\t{ins_nt_dict[Pos][insertion]}"
                                )
                                ntcall_lines["line"][insert_position] += (
                                    f"\t{iabund:.3f}"
                                )
                                if "fs" in i_AAs:
                                    ntcall_lines["line"][insert_position] += "\t\tfs"
                                else:
                                    split_AAs = i_AAs.split(str(AA_pos))
                                    if split_AAs[0]:
                                        ntcall_lines["line"][insert_position] += (
                                            f"\t\t{split_AAs[0]}->{split_AAs[1][:-1]}"
                                        )
                                    else:
                                        ntcall_lines["line"][insert_position] += (
                                            f"\t\t{split_AAs[1]}"
                                        )
                                if (ins_nt_dict[Pos][insertion] >= args.min_count) and (
                                    iabund >= args.min_samp_abund
                                ):
                                    ntcall_lines["variant"][insert_position] = 1
                                i += 1

                        Pos_calls = {}
                        for key in nt_call_dict_dict[Pos]:
                            Pos_calls[key] = nt_call_dict_dict[Pos][key]
                        sorted_calls = sorted(
                            Pos_calls,
                            key=Pos_calls.__getitem__,
                            reverse=True,
                        )

                        ntcall_lines["line"][Pos] = (
                            str(Pos)
                            + "\t"
                            + ref[1][Pos - 1]
                            + "\t"
                            + str(((Pos - 1) // 3) + 1)
                            + "\t"
                            + ref[3][1][((Pos - 1) // 3)]
                        )
                        ntcall_lines["line"][Pos] += (
                            "\t"
                            + str(nt_call_dict_dict[Pos]["A"])
                            + "\t"
                            + str(nt_call_dict_dict[Pos]["T"])
                            + "\t"
                            + str(nt_call_dict_dict[Pos]["C"])
                        )
                        ntcall_lines["line"][Pos] += (
                            "\t"
                            + str(nt_call_dict_dict[Pos]["G"])
                            + "\t"
                            + str(nt_call_dict_dict[Pos]["-"])
                            + "\t"
                            + str(nt_call_dict_dict[Pos]["N"])
                        )
                        ntcall_lines["line"][Pos] += (
                            "\t"
                            + str(total)
                            + "\t"
                            + sorted_calls[0]
                            + "\t"
                            + str(nt_call_dict_dict[Pos][sorted_calls[0]])
                        )
                        ntcall_lines["line"][Pos] += (
                            f"\t{(nt_call_dict_dict[Pos][sorted_calls[0]]/total):.3f}"
                        )

                        consensus[Pos] = sorted_calls

                for Pos in sorted_Pos:
                    try:
                        ntcall_lines["line"][Pos]
                    except:
                        pass
                    else:
                        if consensus[Pos][0] != ref[1][Pos - 1]:
                            ntcall_lines["variant"][Pos] = 1
                            mod = (Pos) % 3

                            if mod == 0:
                                try:
                                    codon = (
                                        consensus[Pos - 2][0]
                                        + consensus[Pos - 1][0]
                                        + consensus[Pos][0]
                                    )
                                except:
                                    codon = "NNN"
                            elif mod == 2:
                                try:
                                    codon = (
                                        consensus[Pos - 1][0]
                                        + consensus[Pos][0]
                                        + consensus[Pos + 1][0]
                                    )
                                except:
                                    codon = "NNN"
                            elif mod == 1:
                                try:
                                    codon = (
                                        consensus[Pos][0]
                                        + consensus[Pos + 1][0]
                                        + consensus[Pos + 2][0]
                                    )
                                except:
                                    codon = "NNN"
                            ntcall_lines["line"][Pos] += (
                                "\t"
                                + aa_call(codon)
                                + "\t"
                                + singlet_codon_call(Pos, consensus[Pos][0], ref[1])[1]
                            )
                        else:
                            ntcall_lines["line"][Pos] += "\t\t"

                        if (
                            nt_call_dict_dict[Pos][consensus[Pos][1]] >= args.min_count
                        ) and (
                            (nt_call_dict_dict[Pos][consensus[Pos][1]] / total)
                            >= args.min_samp_abund
                        ):
                            ntcall_lines["variant"][Pos] = 1
                            ntcall_lines["line"][Pos] += (
                                f"\t{consensus[Pos][1]}\t{nt_call_dict_dict[Pos][consensus[Pos][1]]}\t{(nt_call_dict_dict[Pos][consensus[Pos][1]]/total):.3f}"
                                + "\t"
                                + singlet_codon_call(Pos, consensus[Pos][1], ref[1])[1]
                            )

                            if (
                                nt_call_dict_dict[Pos][consensus[Pos][2]]
                                >= args.min_count
                            ) and (
                                nt_call_dict_dict[Pos][consensus[Pos][2]] / total
                                >= args.min_samp_abund
                            ):
                                ntcall_lines["variant"][Pos] = 1
                                ntcall_lines["line"][Pos] += (
                                    f"\t{consensus[Pos][2]}\t{nt_call_dict_dict[Pos][consensus[Pos][2]]}\t{(nt_call_dict_dict[Pos][consensus[Pos][2]]/total):.3f}\t{singlet_codon_call(Pos, consensus[Pos][2], ref[1])[1]}"
                                )

                for Pos in ntcall_lines["line"]:
                    ntcall_fh.write(ntcall_lines["line"][Pos])
                    ntcall_fh.write("\n")
                    if args.ntvar == 1:
                        try:
                            ntcall_lines["variant"][Pos]
                        except:
                            pass
                        else:
                            ntcallv_fh.write(ntcall_lines["line"][Pos])
                            ntcallv_fh.write("\n")
                if args.ntvar == 1:
                    ntcallv_fh.close()

            else:
                ntcall_fh.write(
                    "Position\tref NT\tA\tT\tC\tG\t-\tN\t\tTotal\tPrimary NT\tCounts\tAbundance\tSecondary NT\tCounts\tAbundance\tTertiary NT\tCounts\tAbundance\n",
                )
                if args.ntvar == 1:
                    ntcallv_fh.write(
                        "Position\tref NT\tA\tT\tC\tG\t-\tN\t\tTotal\tPrimary NT\tCounts\tAbundance\tSecondary NT\tCounts\tAbundance\tTertiary NT\tCounts\tAbundance\n",
                    )

                for Pos in sorted_Pos:
                    try:
                        total = coverage[Pos]
                    except:
                        total = 0
                        print(f"Coverage of position {Pos} not found")

                    if (
                        total >= (sam_read_count * args.ntabund) or args.wgs == 1
                    ) and total >= args.ntcover:
                        try:
                            ins_nt_dict[Pos]
                        except:
                            pass
                        else:
                            i = 1
                            for insertion in ins_nt_dict[Pos]:
                                insert_position = Pos + (i / 1000)
                                i_nts = insertion.split("insert")[1]
                                AA_pos = ((Pos - 1) // 3) + 1
                                iabund = ins_nt_dict[Pos][insertion] / total
                                ntcall_lines["line"][insert_position] = (
                                    f"{Pos}\t-\t\t\t\t\t\t\t{total}\t{i_nts}\t{ins_nt_dict[Pos][insertion]}"
                                )
                                ntcall_lines["line"][insert_position] += (
                                    f"\t{iabund:.3f}"
                                )

                                if (ins_nt_dict[Pos][insertion] >= args.min_count) and (
                                    iabund >= args.min_samp_abund
                                ):
                                    ntcall_lines["variant"][insert_position] = 1
                                i += 1

                        Pos_calls = {}
                        for key in nt_call_dict_dict[Pos]:
                            Pos_calls[key] = nt_call_dict_dict[Pos][key]
                        sorted_calls = sorted(
                            Pos_calls,
                            key=Pos_calls.__getitem__,
                            reverse=True,
                        )

                        ntcall_lines["line"][Pos] = str(Pos) + "\t" + ref[1][Pos - 1]
                        ntcall_lines["line"][Pos] += (
                            "\t"
                            + str(nt_call_dict_dict[Pos]["A"])
                            + "\t"
                            + str(nt_call_dict_dict[Pos]["T"])
                            + "\t"
                            + str(nt_call_dict_dict[Pos]["C"])
                            + "\t"
                            + str(nt_call_dict_dict[Pos]["G"])
                            + "\t"
                            + str(nt_call_dict_dict[Pos]["-"])
                            + "\t"
                            + str(nt_call_dict_dict[Pos]["N"])
                        )
                        ntcall_lines["line"][Pos] += (
                            "\t"
                            + str(total)
                            + "\t"
                            + sorted_calls[0]
                            + "\t"
                            + str(nt_call_dict_dict[Pos][sorted_calls[0]])
                            + "\t"
                            + f"{(nt_call_dict_dict[Pos][sorted_calls[0]]/total):.3f}"
                        )
                        if ref[1][Pos - 1] != sorted_calls[0]:
                            ntcall_lines["variant"][Pos] = 1
                        if (
                            nt_call_dict_dict[Pos][sorted_calls[1]] >= args.min_count
                        ) and (
                            nt_call_dict_dict[Pos][sorted_calls[1]] / total
                        ) >= args.min_samp_abund:
                            ntcall_lines["variant"][Pos] = 1
                            ntcall_lines["line"][Pos] += (
                                "\t"
                                + sorted_calls[1]
                                + "\t"
                                + str(nt_call_dict_dict[Pos][sorted_calls[1]])
                                + "\t"
                                + f"{(nt_call_dict_dict[Pos][sorted_calls[1]]/total):.3f}"
                            )
                            if (
                                nt_call_dict_dict[Pos][sorted_calls[2]] > args.min_count
                            ) and (
                                nt_call_dict_dict[Pos][sorted_calls[2]] / total
                                > args.min_samp_abund
                            ):
                                ntcall_lines["variant"][Pos] = 1
                                ntcall_lines["line"][Pos] += (
                                    "\t"
                                    + sorted_calls[2]
                                    + "\t"
                                    + str(nt_call_dict_dict[Pos][sorted_calls[2]])
                                    + "\t"
                                    + f"{(nt_call_dict_dict[Pos][sorted_calls[2]]/total):.3f}"
                                )

                for Pos in ntcall_lines["line"]:
                    ntcall_fh.write(ntcall_lines["line"][Pos])
                    ntcall_fh.write("\n")
                    if args.ntvar == 1:
                        try:
                            ntcall_lines["variant"][Pos]
                        except:
                            pass
                        else:
                            ntcallv_fh.write(ntcall_lines["line"][Pos])
                            ntcallv_fh.write("\n")

            ntcall_fh.close()
            if args.ntvar == 1:
                ntcallv_fh.close()

    print(f"End {file} main output")
