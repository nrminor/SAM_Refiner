def covar_deconv(args, samp, covariance_dict, sequence_dict):
    """
    Called to perform chimera removal based on the covars output and print the covar deconv output
    Parameters:
    args - argument values
    samp - sam sample name
    covariance_dict - dict of covar lines
    sequence_dict - dict of unique seq lines

    Functionality:
    A sequence is first determined to likely be a true or chimeric sequence.  The ratio of the frequency of a given covariant sequence
    relative to the product of the abundances of each individual polymorphism that is present in that covariant sequence is calculated.
    If the ratio of the sequence abundance to the product is equal to or greater than 1 (beta), that covariant passes the check.
    Any sequence that has an abundance of 0.3 or greater is automatically passed.  Then the passed sequences, in order of greatest
    observed/expected ratio to least, are assigned a new occurrence count based on their constituent individual polymorphisms. The count
    of the least abundant individual polymorphism is assigned to the sequence and constituent polymorphisms making up the sequence have
    their count reduced by the amount of the least abundant polymorphism. Any sequence not yet processed in which that polymorphism is
    present is removed. This process is repeated until all sequences have been reassessed or removed.

    Returns nothing
    """

    passedseqs = {}
    preservedseqs = {}
    covar_dict = covariance_dict
    for seq in sequence_dict:  # pass check for actual : expected abundance
        if seq != "total" and seq != "singles":
            splitseq = seq.split(" ")
            abund = 1
            for sing in splitseq:
                try:
                    abund = abund * (covar_dict[sing] / covar_dict["total"])
                except:
                    abund = abund * (sequence_dict[seq] / sequence_dict["total"])

            try:
                covarabund = covar_dict[seq] / covar_dict["total"]
            except:
                covarabund = sequence_dict[seq] / sequence_dict["total"]
                covar_dict[seq] = sequence_dict[seq]

            if covarabund >= args.autopass:
                preservedseqs[seq] = max(1, args.beta, (covarabund / abund))
            elif covarabund >= abund * args.beta:
                passedseqs[seq] = covarabund / abund
            elif len(seq.split(" ")) == 1:
                passedseqs[seq] = max(1, args.beta)

    if args.min_samp_abund < 1:
        min_count = args.min_samp_abund * covar_dict["total"]
    else:
        min_count = args.min_samp_abund

    if args.pass_out == 1:  # write passed covars to file if enabled
        fh_pass = open(samp + "_covar_pass.tsv", "w")
        fh_pass.write(
            f"{samp}({covar_dict['total']})\nSequences\tCount\tAbundance\tPass Ratio\n"
        )
        for seq in preservedseqs:
            if covar_dict[seq] >= min_count:
                fh_pass.write(
                    f"{seq}\t{covar_dict[seq]}\t{(covar_dict[seq]/covar_dict['total']):.3f}\t{preservedseqs[seq]}*\n"
                )
        for seq in passedseqs:
            if covar_dict[seq] >= min_count:
                fh_pass.write(
                    f"{seq}\t{covar_dict[seq]}\t{(covar_dict[seq]/covar_dict['total']):.3f}\t{passedseqs[seq]}\n"
                )
        fh_pass.close()

    # sort passed covars
    lensortedpassed = sorted(
        passedseqs, key=lambda key: len(key.split(" ")), reverse=True
    )
    ratiolensortedpassed = sorted(
        lensortedpassed, key=lambda key: passedseqs[key], reverse=True
    )
    sortedsingles = sorted(covar_dict["singles"], key=covar_dict["singles"].__getitem__)
    deconved = {}
    for seq in ratiolensortedpassed:  # reassign counts
        singles = seq.split(" ")
        first = 0
        rem_count = 0
        for sing in sortedsingles:
            if sing in singles:
                if covar_dict["singles"][sing] > 0:
                    if first == 0:
                        first = 1
                        rem_count = covar_dict["singles"][sing]
                        covar_dict["singles"][sing] = 0
                        deconved[seq] = rem_count
                    else:
                        covar_dict["singles"][sing] = (
                            covar_dict["singles"][sing] - rem_count
                        )
                else:
                    break
        sortedsingles = sorted(
            covar_dict["singles"], key=covar_dict["singles"].__getitem__
        )

    sortedpreserved = sorted(preservedseqs, key=lambda key: covar_dict[key])

    for seq in sortedpreserved:
        singles = seq.split(" ")
        first = 0
        rem_count = 0
        for sing in sortedsingles:
            if sing in singles:
                if covar_dict["singles"][sing] > 0:
                    if first == 0:
                        first = 1
                        rem_count = covar_dict["singles"][sing]
                        covar_dict["singles"][sing] = 0
                        deconved[seq] = rem_count
                    else:
                        covar_dict["singles"][sing] = (
                            covar_dict["singles"][sing] - rem_count
                        )
                else:
                    break
        sortedsingles = sorted(
            covar_dict["singles"], key=covar_dict["singles"].__getitem__
        )

    newtotal = sum(deconved.values())
    fh_deconv = open(samp + "_covar_deconv.tsv", "w")
    fh_deconv.write(
        f"{samp}({covar_dict['total']}) | ({newtotal})\nSequences\tCount\tAbundance\n"
    )
    sorted_deconved = sorted(deconved, key=deconved.__getitem__, reverse=True)
    for seq in sorted_deconved:  # write deconv
        if deconved[seq] >= min_count:
            fh_deconv.write(f"{seq}\t{deconved[seq]}\t{(deconved[seq]/newtotal):.3f}\n")
    fh_deconv.close()

    print(f"End covar deconv out for {samp}")  # END COVAR DECONV OUT

    return ()


def dechim(args, seqs):
    """
    Called to perform chimera removal based on the unique seqs output
    Parameters:
    args - argument values
    seqs - dict of unique seq lines
    Functionality:
    The individual unique sequences, starting with the lowest abundance, are broken up into all possible dimeric halves. Each pair
    is then compared to all the other sequences to detect potential parents. A sequence is flagged as a potential parent if its
    abundance is greater than or equal to the abundance of the potential chimera multiplied by foldab and there is at least one
    other sequence that would be a matched parent to the complimentary dimeric half. When a pair of dimeric halves have potential
    parents, the abundances of parent pairs are multiplied. The products from each potential parent pairings are summed as an
    expected abundance value and compared to the observed abundance of the potential chimera. If the abundance of the potential
    chimera is less than that of the expected value multiplied by alpha, that sequence is flagged as a chimera and removed. The
    counts attributed to the flagged chimeric sequence are then redistributed to the parent sequences based on the relative expected
    contribution to recombination.
    Returns new sequence dict without found chimeras
    """

    total = seqs["total"]
    del seqs["total"]
    sorted_seqs = sorted(
        seqs, key=seqs.__getitem__
    )  # sort sequences by abundance, least to greatest
    chimeras = []
    for seq in sorted_seqs:
        pot_chim = ["Reference"] + seq.split() + ["Reference"]
        chim_halves = []
        for i in range(len(pot_chim) - 1):  # create dimera halves
            chim_halves.append([pot_chim[: i + 1], pot_chim[i + 1 :]])
        parent_pairs = []

        for dimera in chim_halves:  # check for potential parents
            pot_left = []
            pot_rt = []
            lft_len = len(dimera[0])
            rt_len = len(dimera[1])
            for pseq in sorted_seqs:
                if seq != pseq:
                    if seqs[pseq] >= (seqs[seq] * args.foldab):
                        pot_par = ["Reference"] + pseq.split() + ["Reference"]
                        if dimera[0] == pot_par[:lft_len]:
                            pot_left.append(pot_par)
                        if (len(pot_par) >= rt_len) and (
                            dimera[1] == pot_par[(len(pot_par) - rt_len) :]
                        ):
                            pot_rt.append(pot_par)

            if len(pot_left) > 0 and len(pot_rt) > 0:
                for left_par in pot_left:  # check potential parents' pairing
                    for rt_par in pot_rt:
                        if left_par != rt_par:
                            left_break = left_par[lft_len]
                            rt_break = rt_par[(len(rt_par) - rt_len) - 1]
                            if left_break == "Reference" or rt_break == "Reference":
                                parent_pairs.append(
                                    [" ".join(left_par[1:-1]), " ".join(rt_par[1:-1])]
                                )
                            else:
                                left_break_Pos = ""
                                for c in left_break:
                                    if c.isdigit():
                                        left_break_Pos += c
                                    elif left_break_Pos:
                                        break

                                rt_break_Pos = ""
                                for c in rt_break:
                                    if c.isdigit():
                                        rt_break_Pos += c
                                    elif rt_break_Pos:
                                        break

                                if int(left_break_Pos) > int(rt_break_Pos):
                                    parent_pairs.append(
                                        [
                                            " ".join(left_par[1:-1]),
                                            " ".join(rt_par[1:-1]),
                                        ]
                                    )

        par_tot_abund = 0
        pair_probs = []
        for parents in parent_pairs:  # calc 'expected' abundance
            pair_prob = (seqs[parents[0]] / total) * (seqs[parents[1]] / total)
            par_tot_abund += pair_prob
            pair_probs.append(pair_prob)

        recomb_count = par_tot_abund * total

        if not seqs[seq] >= recomb_count * args.alpha:  # chimera check
            redist_count = float(seqs[seq])
            chimeras.append(seq)
            seqs[seq] = 0
            if args.redist == 1:  # redist counts of chimera
                toadd = {}
                for i in range(len(parent_pairs)):
                    counts_to_redist = (
                        redist_count * (pair_probs[i] / par_tot_abund)
                    ) / 2
                    seqs[parent_pairs[i][0]] += counts_to_redist
                    seqs[parent_pairs[i][1]] += counts_to_redist

    for chim in chimeras:  # remove chimeras
        del seqs[chim]

    total = sum(seqs.values())

    seqs["total"] = total

    return seqs


def chim_rm(args, samp, sequences):
    """
    Called to reiterate the chim_rm chimera removal process based on unique seqs and print the chim_rm output
    Parameters:
    args - argument values
    samp - sample name from sam file
    sequencess - dict of unique sequences
    Functionality:
    Looped calling of chimera removal until no more chimeras are found or limit hit.  Results are then printed
    Returns nothing
    """

    seqs = sequences
    pre_len = len(seqs)
    inf_loop_shield = 0
    while True:  # send sequences for chimera removal while chimeras are still found
        dechim(args, seqs)
        post_len = len(seqs)
        inf_loop_shield += 1
        if post_len >= pre_len:
            break
        if inf_loop_shield > args.max_cycles:
            break
        pre_len = len(seqs)

    total = seqs["total"]
    del seqs["total"]
    if args.min_samp_abund < 1:
        min_count = args.min_samp_abund * total
    else:
        min_count = args.min_samp_abund
    sorted_seqs = sorted(seqs, key=seqs.__getitem__, reverse=True)
    fh_dechime = open(
        f"{samp}_a{args.alpha}f{args.foldab}rd{args.redist}_chim_rm.tsv", "w"
    )
    fh_dechime.write(f"{samp}({int(total)})\n")
    fh_dechime.write("Sequences\tCount\tAbundance\n")
    for seq in seqs:  # write chim_rm seqs
        abund = seqs[seq] / total
        if seqs[seq] >= min_count:
            fh_dechime.write(f"{seq}\t{round(seqs[seq])}\t{abund:.3f}\n")

    fh_dechime.close()
    print(f"End chim_rm out for {samp}")
    return ()


def chim_process(args, samp):
    """
    Called to collect sequences from covar and unique seq files, then pass those to the chimera removal functions
    Parameters:
    args - argument values
    samp - samples name from sam file
    Functionality:
    Tries to open files based on based sample name.  If successful, calls the specific chimera removal
    functions
    Returns nothing
    """
    in_seqs = {}
    try:
        seqin_file = open(samp + "_unique_seqs.tsv")
        for line in seqin_file:
            lineparts = line.strip("\n\r").split("\t")
            try:
                lineparts[1]
            except:
                in_seqs = {"total": int(lineparts[0].split("(")[1][0:-1])}
            else:
                if lineparts[1] != "Count":
                    if float(lineparts[2]) >= args.chim_in_abund:
                        in_seqs[lineparts[0]] = float(lineparts[1])
        seqin_file.close()
    except:
        print(f"Reading of {samp}_unique_seqs.tsv failed")

    if args.deconv == 1:
        in_covars = {}
        try:
            covin_file = open(samp + "_covars.tsv")
            for line in covin_file:
                lineparts = line.strip("\n\r").split("\t")
                try:
                    lineparts[1]
                except:
                    in_covars = {
                        "total": int(lineparts[0].split("(")[1][0:-1]),
                        "singles": {},
                    }
                else:
                    if lineparts[1] != "Count":
                        if float(lineparts[2]) >= args.chim_in_abund:
                            in_covars[lineparts[0]] = int(lineparts[1])
                            if len(lineparts[0].split(" ")) == 1:
                                in_covars["singles"][lineparts[0]] = int(lineparts[1])
            covin_file.close()
            if in_covars and in_seqs:
                covar_deconv(args, samp, in_covars, in_seqs)
        except:
            print(f"Reading of {samp}_covars.tsv failed")

    if args.chim_rm == 1:
        chim_rm(args, samp, in_seqs)
