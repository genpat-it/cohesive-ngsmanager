#!/usr/bin/env python3

import os
import sys

def analyse_innuendo(
        # sample, dictss, pathr,
        summary_file,
        sampleds,
        cmp_cd,
        base,
        virulence_dict, resfinder_dict, card_dict, ncbi_dict, argannot_dict,
        outputfileamr, outputfilevf, outputdedotto, listclasses):
    res_innu = open(summary_file, "r").readlines()

    # cmp_cd = dictss[sampleds].replace("-", ".")
    headeratb = res_innu[0].strip().split("\t")
    # for row in res_innu:
    dizatb = {}
    for classe in listclasses:
        dizatb[classe] = []
    dizunknown = []
    # if sample in row:
    row = res_innu[1]
    indiceid = 1
    listgenes = (row.strip().split("\t")[2:])
    for geneid in listgenes:
        indiceid += 1
        if geneid != ".":
            if ";" in geneid:
                gl = []
                tmpminilist = geneid.split(";")
                for g in tmpminilist:
                    gl.append(float(g))
                geneid = max(gl)
            if float(geneid) >= 0:  # da modificare se vogliamo mettere thresholds
                geneori = headeratb[indiceid]
                gene = geneori.lower().replace("escherichia_coli_", "escherichia coli ").replace("(mls)lin",
                                                                                                 "lin").split(
                    "_")[0].split(".")[0]
                if "escherichia" in gene or "listeria" in gene or "brucella" in gene or "campylobacter" in gene:
                    for chiave in card_dict.keys():
                        if gene in chiave:
                            if card_dict[chiave][0] not in dizatb.keys():
                                dizatb[card_dict[chiave][0]] = [geneori]
                            else:
                                dizatb[card_dict[chiave][0]].append(geneori)
                elif gene.split("(")[0] in virulence_dict.keys() or gene in virulence_dict.keys():
                    if gene.split("(")[0] in virulence_dict.keys():
                        if virulence_dict[gene.split("(")[0]] not in dizatb.keys():
                            dizatb[virulence_dict[gene.split("(")[0]]] = [geneori]
                        else:
                            dizatb[virulence_dict[gene.split("(")[0]]].append(geneori)
                    else:
                        if virulence_dict[gene][0] not in dizatb.keys():
                            dizatb[virulence_dict[gene]] = [geneori]
                        else:
                            dizatb[virulence_dict[gene]].append(geneori)
                elif gene.split("(")[0] in card_dict.keys() or gene in card_dict.keys():
                    if gene.split("(")[0] in card_dict.keys():
                        if card_dict[gene.split("(")[0]][0] not in dizatb.keys():
                            dizatb[card_dict[gene.split("(")[0]][0]] = [geneori]
                        else:
                            dizatb[card_dict[gene.split("(")[0]][0]].append(geneori)
                    else:
                        if card_dict[gene][0] not in dizatb.keys():
                            dizatb[card_dict[gene][0]] = [geneori]
                        else:
                            dizatb[card_dict[gene][0]].append(geneori)
                elif gene.split("(")[0] in resfinder_dict.keys() or gene in resfinder_dict.keys():
                    if gene.split("(")[0] in resfinder_dict.keys():
                        if resfinder_dict[gene.split("(")[0]][0] not in dizatb.keys():
                            dizatb[resfinder_dict[gene.split("(")[0]][0]] = [geneori]
                        else:
                            dizatb[resfinder_dict[gene.split("(")[0]][0]].append(geneori)
                    else:
                        if resfinder_dict[gene][0] not in dizatb.keys():
                            dizatb[resfinder_dict[gene][0]] = [geneori]
                        else:
                            dizatb[resfinder_dict[gene][0]].append(geneori)
                elif gene in argannot_dict.keys():
                    if argannot_dict[gene][0] not in dizatb.keys():
                        dizatb[argannot_dict[gene][0]] = [geneori]
                    else:
                        dizatb[argannot_dict[gene][0]].append(geneori)
                elif gene in ncbi_dict.keys():
                    if ncbi_dict[gene][0] not in dizatb.keys():
                        dizatb[ncbi_dict[gene][0]] = [geneori]
                    else:
                        dizatb[ncbi_dict[gene][0]].append(geneori)
                else:
                    dizunknown.append(geneori)
    outputfileamr.write(cmp_cd + "," + sampleds.replace("DS", "") + ",")
    outputfilevf.write(cmp_cd + "," + sampleds.replace("DS", "") + ",")
    outputdedotto.write(cmp_cd + "," + sampleds.replace("DS", "") + ",")
    almeno1file = open(base + "_abricate_AMR.check", "w")
    almeno1filevf = open(base + "_abricate_VF.check", "w")
    almeno1 = 0
    almeno1vf = 0
    list_amr = []
    for classe in listclasses:
        if classe != "virulence" and classe != "efflux" and classe != "others" and classe != "unknown":
            outputfileamr.write(";".join(sorted(set(dizatb[classe]))) + ",")
            if len(set(dizatb[classe])) >= 1:
                list_amr.append(classe)
                almeno1 = 1
        if classe == "virulence":
            outputfilevf.write(";".join(dizatb[classe]) + ",")
            if len(set(dizatb[classe])) >= 1:
                almeno1vf = 1
        if classe == "efflux":
            outputfileamr.write(";".join(dizatb[classe]) + ",")
            if len(set(dizatb[classe])) >= 1:
                list_amr.append(classe)
                almeno1 = 1
        if classe == "others":
            outputfileamr.write(";".join(dizatb[classe]) + ",")
            if len(set(dizatb[classe])) >= 1:
                list_amr.append(classe)
                almeno1 = 1
        if classe == "unknown":
            outputfileamr.write(";".join(dizunknown) + ";".join(dizatb[classe]) + ",")
            if len(set(dizatb[classe])) >= 1:
                list_amr.append(classe)
                almeno1 = 1
    outputdedotto.write(";".join(list_amr))
    outputdedotto.write(",")
    if almeno1vf == 1:
        outputdedotto.write("virulence,")
    outputfileamr.write("\n")
    outputfilevf.write("\n")
    if almeno1 == 0:
        almeno1file.write("FAIL")
        outputdedotto.write("FAIL,")
    if almeno1 == 1:
        almeno1file.write("PASS")
        outputdedotto.write("PASS,")
    if almeno1vf == 0:
        almeno1filevf.write("FAIL")
        outputdedotto.write("FAIL")
    if almeno1vf == 1:
        almeno1filevf.write("PASS")
        outputdedotto.write("PASS")
    almeno1file.close()
    almeno1filevf.close()
    outputdedotto.write("\n")


def dict_resfinder():
    resfinder_dict = {}
    listfolder = os.listdir("/mnt/biowork/databases/PROGRAMS/resfinder/resfinder_db/")
    for fileatb in listfolder:
        if ".fsa" in fileatb:
            if "macrolide" in fileatb:
                atb = "MLS"
            elif "sulphonamide" in fileatb:
                atb = "trimethoprim"
            elif "beta-lactam" in fileatb:
                atb = "betalactam"
            else:
                atb = fileatb.split(".fsa")[0]
            atb_genes = open("/mnt/biowork/databases/PROGRAMS/resfinder/resfinder_db/" + fileatb).readlines()
            for riga in atb_genes:
                if riga[0] == ">":
                    gene = riga.split("_")[0][1:]
                    if gene not in resfinder_dict.keys():
                        resfinder_dict[gene.strip().lower()] = []
                        resfinder_dict[gene.strip().lower()].append(atb)
                    else:
                        resfinder_dict[gene.strip().lower()].append(atb)
    listfolderpf = os.listdir("/mnt/biowork/databases/PROGRAMS/plasmidfinder/plasmidfinder_db/")
    for fileclass in listfolderpf:
        if ".fsa" in fileclass:
            pl_genes = open("/mnt/biowork/databases/PROGRAMS/plasmidfinder/plasmidfinder_db/" + fileclass).readlines()
            for rigapl in pl_genes:
                if rigapl[0] == ">":
                    gene = rigapl.split("_")[0][1:]
                    resfinder_dict[gene.lower()] = "others"
    return resfinder_dict


def dict_ncbi(diziotranslate):
    ncbi_dict = {}
    atb_genes = open("/mnt/biowork/databases/PROGRAMS/abricate/ncbi/sequences", "r").readlines()
    for riga in atb_genes:
        if riga[0] == ">":
            gene = riga.split("~~~")[1].lower()
            atb = diziotranslate[riga.split("~~~")[-1].split(" ")[0].split("/")[0].lower()]
            if gene not in ncbi_dict.keys():
                ncbi_dict[gene.strip().lower()] = []
                ncbi_dict[gene.strip().lower()].append(atb)
            else:
                ncbi_dict[gene.strip().lower()].append(atb)
    return ncbi_dict


def dict_argannot(diziotranslate):
    argannot_dict = {}
    atb_genes = open("/mnt/biowork/databases/PROGRAMS/abricate/argannot/sequences", "r").readlines()
    for riga in atb_genes:
        if riga[0] == ">":
            gene = riga.split("~~~")[1].lower()
            atb = diziotranslate[riga.split("~~~")[-1].split("(")[1].split(")")[0].split("_")[0].lower()]
            if gene not in argannot_dict.keys():
                argannot_dict[gene.strip().lower()] = []
                argannot_dict[gene.strip().lower()].append(atb)
            else:
                argannot_dict[gene.strip().lower()].append(atb)
    return argannot_dict


def create_diziotranslate():
    danielfile = open("/mnt/biowork/databases/PROGRAMS/abricate/ListofAB.txt").readlines()
    diziotranslate = {}
    for riga in danielfile:
        classe = riga.split(":")[0].strip()
        atbs = riga.split(":")[1].strip().split(",")
        for atb in atbs:
            atb = atb.split()[0].lower()
            if atb not in diziotranslate.keys():
                diziotranslate[atb.strip()] = classe
    return diziotranslate


def dict_card(diziotranslate):
    obo = open("/mnt/biowork/databases/PROGRAMS/card/aro.obo").readlines()
    name = ""
    card_dict = {}
    i = 1
    lista = []
    for _ in [1, 2]:
        for riga in obo:
            if riga.strip() == "[Term]":
                i = 0
                lista = []
            if "name:" in riga:
                name = riga.split(":")[1].split("with")[0].strip().lower()
                if "antibiotic efflux pump" in name:
                    card_dict[name] = ["efflux"]
            if "relationship: confers_resistance_to" in riga:
                drug = riga.split("!")[1].split()[0].strip().lower()
                if drug in diziotranslate.keys():
                    if diziotranslate[drug] not in lista:
                        lista.append(diziotranslate[drug])
                i = 1
            if "is_a: " in riga:
                name2 = riga.split("!")[1].split("with")[0].strip().lower()
                if "efflux" in riga:
                    card_dict[name] = ["efflux"]
                if name2 in card_dict.keys():
                    card_dict[name] = card_dict[name2]
                    i = 0
            if i == 1 and len(lista) > 0:
                card_dict[name] = lista
        for ele in card_dict.keys():
            if len(card_dict[ele]) > 1:
                if "colistin" in card_dict[ele] and 'polypeptides' in card_dict[ele]:
                    card_dict[ele].remove('polypeptides')
                else:
                    card_dict[ele] = ["efflux"]
    return card_dict


def dict_vfdb():
    virulence_dict = {}
    list_vfdb = open("/mnt/biowork/databases/PROGRAMS/vfdb/lista_VFDB").readlines()
    for riga in list_vfdb:
        virulence_dict[riga.strip().lower()] = "virulence"
    listfoldervf = os.listdir("/mnt/biowork/databases/PROGRAMS/virulencefinder/virulencefinder_db/")
    for filepathogen in listfoldervf:
        if ".fsa" in filepathogen:
            vl_genes = open(
                "/mnt/biowork/databases/PROGRAMS/virulencefinder/virulencefinder_db/" + filepathogen).readlines()
            for rigavl in vl_genes:
                if rigavl[0] == ">":
                    genev = rigavl.split("_")[0].split(":")[0][1:]
                    virulence_dict[genev.lower()] = "virulence"
    return virulence_dict


def main():
    if len(sys.argv) < 5:
        sys.exit(1)

    summary_file = sys.argv[1]
    ds = sys.argv[2]
    cmp = sys.argv[3]
    dt = sys.argv[4]

    base = f"{ds}-{dt}_{cmp}"

    diziotranslate = create_diziotranslate()
    card_dict = dict_card(diziotranslate)
    ncbi_dict = dict_ncbi(diziotranslate)
    argannot_dict = dict_argannot(diziotranslate)
    resfinder_dict = dict_resfinder()
    virulence_dict = dict_vfdb()
    listclasses = ["aminoglycoside", "betalactam", "polypeptides", "colistin", "quinolone", "fosfomycin",
                   "fusidicacid", "glycopeptide", "MLS", "nitroimidazole", "oxazolidinone", "phenicol",
                   "rifampicin", "trimethoprim", "tetracycline", "others", "efflux", "unknown", "virulence"]
    # path = "/bioinfonas/analisi/run_results/"
    errori = open(base + "_errorianalisi.txt", "w")
    # run = sys.argv[1]

    outputdedotto = open(base + "_import_abricate.csv", "w")
    outputdedotto.write("CMP_CD,SAMPLE_DS,AMR,VIRULF,CHECK_AMR,CHECK_VIRULF\n")
    # dictss = {}
    # check = 0
    # SS = open(path + run + "/SampleSheet.csv", "r").readlines()
    # for rigass in SS:
    #     if "sample_id" in rigass.lower() and check == 0:
    #         check = 1
    #     elif check == 1:
    #         dictss[rigass.split(",")[1].split("-")[0]] = rigass.split(",")[0]
    # for cartella in os.listdir(path + run):
    #     if os.path.isdir(path + run + "/" + cartella) == True:
    #         pathr = path + run + "/" + cartella + "/"
    #         for elemento in os.listdir(pathr):
    #             if "_spades_scaffolds_L200.fasta" in elemento:
    #                 key = elemento.replace("_spades_scaffolds_L200.fasta", "")
    outputfileamr = open(base + "_output_abricate_AMR.csv", "w")
    outputfilevf = open(base + "_output_abricate_VF.csv", "w")
    outputfileamr.write("CMP_CD,SAMPLE_DS,")
    outputfilevf.write("CMP_CD,SAMPLE_DS,")
    for cla in listclasses:
        if cla == "virulence":
            outputfilevf.write(cla + "\n")
        if cla != "virulence":
            outputfileamr.write(cla + ",")
    outputfileamr.write("\n")
    try:
        analyse_innuendo(summary_file,
                         ds,
                         cmp,
                         base,
                         virulence_dict, resfinder_dict, card_dict, ncbi_dict,
                         argannot_dict, outputfileamr, outputfilevf, outputdedotto, listclasses
                         )
    except Exception as ex:
        errori.write("exception :" + str(ex))
    outputfileamr.close()
    outputfilevf.close()
    outputdedotto.close()


if __name__ == "__main__":
    main()
