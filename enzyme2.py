def reverse(s):
    r = ""
    for c in s:
        r = c + r

    return r

def comp(s):
    r = ""
    for c in s:
        if(c == "A"):
            r = r + "T"
        elif(c == "T"):
            r = r + "A"
        elif(c == "G"):
            r = r + "C"
        elif(c == "C"):
            r = r + "G"
        else:
            r = r + c
    return(r)

code = raw_input("Code: ")

#enzyme = raw_input("Enzyme: ")

#comprev = comp(reverse(enzyme))

#lene = len(enzyme)

#lenc = len(comprev)

start = "\033[1m"

end = "\033[0;0m"

yes = False

no = True

pal = False

listenzyme = [("AACGTT", "AclI"),("AAGCTT", "HindIII"),("AATATT","SspI"),("AATT","MluCI Tsp509I"),("ACATGT","PciI"),("ACCGGT","AgeI AgeI-HF"),("ACCTGC","BspMI BfuAI"),("ACCWGGT","SexAI"),("ACGCGT","MluI"),("ACGGC","BceAI"),("ACGT","HpyCH4IV"),("ACNGT","HpyCH4III"),("ACNNNNGTAYC","BaeI"),("ACNNNNNCTCC", "BsaXI"),("ACRYGT", "AflIII"),("(ACTAGT", "SpeI SpeI-HF"),("ACTGG", "BsrI"),("ACTGGG", "BmrI"),("AGATCT", "BglII"),("AGCGCT ", "AfeI"),("AGCT",  "AluI"),("AGGCCT", "StuI"),("AGTACT", "ScaI ScaI-HF"),("ATCGAT", "ClaI BspDI"),("ATCTATGTCGGGTGCGGAGAAAGAGGTAAT", "PI-SceI"),("ATGCAT", "NsiI NsiI-HF"),("ATTAAT", "AseI"),("ATTTAAAT",   "SwaI"),("CAANNNNNGTGG",  "CspCI"),("CAATTG", "MfeI MfeI-HF"),("CACGAG",   "BssSI BssSaI"),("CACGAG",  "Nb.BssSI"),("CACGTC",   "BmgBI"),("CACGTG", "PmlI"),("CACNNNGTG", "DraIII DraIII-HF"),("CACNNNNGTG", "AleI"),("CAGCAG",  "EcoP15I"),("CAGCTG" ,"PvuII PvuII-HF"),("CAGNNNCTG",  "AlwNI"),("CAGTG",  "BtsIMutI"),("NNCASTGNN",  "TspRI"),("CATATG", "NdeI"),("CATG", "NlaIII"),("CATG",   "CviAII"),("CATG",   "FatI"),("CAYNNNNRTG", "MslI"),("CC",   "FspEI"),("CCANNNNNNNNNTGG","XcmI"),("CCANNNNNNTGG",   "BstXI"),("CCANNNNNTGG",    "PflMI"),("CCATC",  "BccI"),("CCATGG", "NcoI NcoI-HF"),("CCCAGC",   "BseYI"),("CCCGC",  "FauI"),("CCCGGG", "SmaI"),("CCCGGG", "XmaI TspMI"),("CCD",   "Nt.CviPII"),("CCDG", "LpnPI"),("CCGC", "AciI"),("CCGCGG", "SacII"),("CCGCTC",   "BsrBI"),("CCGG" ,  "MspI HpaII"),("CCNGG",  "ScrFI"),("CCNGG" , "BssKI StyD4I"),("CCNNGG", "BsaJI"),("CCNNNNNNNGG",    "BslI"),("CCRYGG", "BtgI"),("CCSGG",  "NciI"),("CCTAGG", "AvrII"),("CCTC" ,  "MnlI"),("CCTCAGC",  "BbvCI"),("CCTCAGC", " Nb.BbvCI"),("CCTCAGC",  "Nt.BbvCI"),("CCTGCAGG",   "SbfI SbfI-HF"),("CCTNAGC",  "Bpu10I"),("CCTNAGG",    "Bsu36I"),("CCTNNNNNAGG",    "EcoNI"),("CCTTC",  "HpyAV"),("CCWGG",  "BstNI"),("CCWGG",  "PspGI"),("CCWWGG", "StyI StyI-HF"),("CGANNNNNNTGC",  "BcgI"),("CGATCG","PvuI PvuI-HF"),("CGCG",   "BstUI"),("CGGCCG", "EagI EagI-HF"),("CGGWCCG",    "RsrII"),("CGRYCG", "BsiEI"),("CGTACG" ,"BsiWI BsiWI-HF"),("CGTCTC" ,"BsmBI"),("CGWCG" , "Hpy99I"),("CMGCKG", "MspA1I"),("CNNR" , "MspJI"),("CRCCGGYG",   "SgrAI"),("CTAG"  , "BfaI"),("CTCAG" , "BspCNI"),("CTCGAG", "XhoI PaeR7I TliI"),("CTCTTC", "EarI"),("CTGAAG",   "AcuI"),("CTGCAG" ,"PstI PstI-HF"),("CTGGAG"  , "BpmI"),("CTNAG" , "DdeI"),("CTRYAG" ,"SfcI"),("CTTAAG", "AflII"),("CTTGAG" ,  "BpuEI"),("CTYRAG", "SmlI"),("CYCGRG" ,"AvaI BsoBI"),("GAAGA"  ,"MboII"),("GAAGAC", "BbsI BbsI-HF"),("GAANNNNTTC", "XmnI"),("GAATGC"   , "BsmI"),("GAATGC" , "Nb.BsmI"),("GAATTC", "EcoRI EcoRI-HF"),("GACGC", "HgaI"),("GACGTC", "AatII"),("GACGTC", "ZraI"),("GACNNNGTC",  "Tth111I PflFI"),("GACNNNNGTC", "PshAI"),("GACNNNNNGTC",    "AhdI"),("GACNNNNNNGTC",   "DrdI"),("GAGCTC", "Eco53kI"),("GAGCTC" ,"SacI SacI-HF"),("GAGGAG" ,   "BseRI"),("GAGTC"  ,"PleI"),("GAGTC" ,"Nt.BstNBI"),("GAGTC"  ,"MlyI"),("GANTC"  ,"HinfI"),("GATATC" ,"EcoRV EcoRV-HF"),("GATC"  , "MboI Sau3AI DpnII BfuCI"),("GATC"  , "DpnI"),("GATNNNNATC", "BsaBI"),("GAWTC" , "TfiI"),("GCAATG" ,"BsrDI"),("GCAATG"  ,"Nb.BsrDI"),("GCAGC", "BbvI"),("GCAGTG", "BtsI BtsaI"),("GCAGTG" , "Nb.BtsI"),("GCANNNNNTGC",    "BstAPI"),("GCATC" , "SfaNI"),("GCATGC", "SphI SphI-HF"),("GCCCGGGC" ,  "SrfI"),("GCCGAG"  , 'NmeAIII'),("GCCGGC" ,"NaeI"),("GCCGGC" ,"NgoMIV"),("GCCNNNNNGGC" ,   "BglI"),("GCGATCGC"  , "AsiSI"),("GCGATG"  , "BtgZI"),("GCGC"  , "HinP1I"),("GCGC"  , "HhaI"),("GCGCGC", "BssHII"),("GCGGCCGC",   "NotI NotI-HF"),("GCNGC" , "Fnu4HI"),("GCNNGC", "Cac8I"),("GCNNNNNNNGC" ,  "MwoI"),("GCTAGC", "NheI NheI-HF"),("GCTAGC" ,"BmtI BmtI-HF"),("GCTCTTC" ,   "SapI BspQI"),("GCTCTTC"  , "Nt.BspQI"),("GCTNAGC"  ,  "BlpI"),("GCWGC" , "TseI ApeKI"),("GDGCHC" ,"Bsp1286I"),("GGATC" , "AlwI"),("GGATC" ,"Nt.AlwI"),("GGATCC", "BamHI BamHI-HF"),("GGATG" ,"FokI"),("GGATG"  ,"BtsCI"),("GGCC"   ,"HaeIII PhoI"),("GGCCGGCC",   "FseI"),("GGCCNNNNNGGCC",  "SfiI"),("GGCGCC", "NarI"),("GGCGCC" ,"KasI"),("GGCGCC" ,"SfoI"),("GGCGCC" ,"PluTI"),("GGCGCGCC" ,  "AscI"),("GGCGGA"   , "EciI"),("GGGAC"   , "BsmFI"),("GGGCCC" ,"ApaI"),("GGGCCC" ,"PspOMI"),("GGNCC" , "Sau96I"),("GGNNCC", "NlaIV"),("GGTACC" ,"KpnI KpnI-HF"),("GGTACC" ,"Acc65I"),("GGTCTC" ,"BsaI BsaI-HF"),("GGTGA"  ,"HphI"),("GGTNACC" ,   "BstEII  BstEII-HF"),("GGWCC" , "AvaII"),("GGYRCC" ,"BanI"),("GKGCMC" ,"BaeGI"),("GRCGYC" ,"BsaHI"),("GRGCYC" ,"BanII"),("GTAC"  , "RsaI"),("GTAC"  , "CviQI"),("GTATAC", "BstZ17I"),("GTATAC" , "BstZ17I-HF"),("GTATCC" ,"BciVI"),("GTCGAC" ,"SalI SalI-HF"),("GTCTC", "Nt.BsmAI"),("GTCTC" , "BsmAI BcoDI"),("GTGCAC", "ApaLI"),("GTGCAG" ,  "BsgI"),("GTMKAC" ,"AccI"),("GTNNAC" ,"Hpy166II"),("GTSAC" , "Tsp45I"),("GTTAAC", "HpaI"),("GTTTAAAC",   "PmeI"),("GTYRAC" ,"HincII"),("GWGCWC", "BsiHKAI"),("RAATTY" ,"ApoI ApoI-HF"),("RCATGY" ,"NspI"),("RCCGGY" ,"BsrFI BsrFaI"),("RGATCY" ,"BstYI"),("RGCGCY" ,"HaeII"),("RGCY"  , "CviKI-1"),("RGGNCCY","EcoO109I"),("RGGWCCY","PpuMI"),("TAACTATAACGGTCCTAAGGTAGCGAA", "I-CeuI"),("TACGTA" ,"SnaBI"),("TAGGGATAACAGGGTAAT" , "I-SceI"),("TCATGA", "BspHI"),("TCCGGA" ,"BspEI"),("TCCRAC"  , "MmeI"),("TCGA"  , "TaqaI"),("TCGCGA", "NruI NruI-HF"),("TCNGA"  ,"Hpy188I"),("TCNNGA", "Hpy188III"),("TCTAGA" ,"XbaI"),("TGATCA" ,"BclI"),("TGCA"  , "HpyCH4V"),("TGCGCA", "FspI"),("TGGCAAACAGCTATTATGGGTATTATGGGT" ,"PI-PspI"),("TGGCCA", "MscI"),("TGTACA" ,"BsrGI BsrGI-HF"),("TTAA"  , "MseI"),("TTAATTAA",   "PacI"),("TTATAA" ,"PsiI"),("TTCGAA" ,"BstBI"),("TTTAAA" ,"DraI"),("VCTCGAGB" ,  "PspXI"),("WCCGGW" ,"BsaWI"),("YACGTR" ,"BsaAI"),("YGGCCR" ,"EaeI")]


for z in listenzyme:
    pal = False
    enzyme = z[0]
    comprev = comp(reverse(enzyme))
    if comprev == enzyme:
        pal = True
    lene = len(enzyme)
    lenc = len(comprev)
    for x in range(0, len(code)-lene):

        for y in range(0, lene):
            if(enzyme[y]=="R" and not ([x+y]=="A" or code[x+y]=="G")):
                    no = False
            elif(enzyme[y]=="Y" and not (code[x+y]=="C" or code[x+y]=="T")):
                 no = False
            elif(enzyme[y]=="M" and not (code[x+y]=="A" or code[x+y]=="C")):    
                no = False
            elif(enzyme[y]=="K" and not (code[x+y]=="G" or code[x+y]=="T")):
                no = False
            elif(enzyme[y]=="S" and not (code[x+y]=="C" or code[x+y]=="G")):
                no = False
            elif(enzyme[y]=="W" and not (code[x+y]=="A" or code[x+y]=="T")):
                no = False
            elif(enzyme[y]=="H" and not (code[x+y]=="A" or code[x+y]=="C" or code[x+y]=="T")):
                no = False
            elif(enzyme[y]=="V" and not (code[x+y]=="A" or code[x+y]=="G" or code[x+y]=="C")):
                no = False
            elif(enzyme[y]=="B" and not (code[x+y]=="C" or code[x+y]=="G" or code[x+y]=="T")):
                no = False
            elif(enzyme[y]=="D" and not (code[x+y]=="A" or code[x+y]=="G" or code[x+y]=="T")):
                no = False
            elif(enzyme[y]=="N" and not (code[x+y]=="A" or code[x+y]=="G" or code[x+y]=="T" or code[x+y]=="C")):
                no = False


            elif (enzyme[y]!=code[x+y]):
                no = False

        if no is True:
            print ''
            print "regular"
            print z[1]
            print enzyme
            print code[:x] + start + code[x:x+lene] + end + code[x+lene:]
        else:
            no = True

        if pal is False:

            for t in range(0, lenc):

                if(comprev[t]=="R" and not ([x+t]=="A" or code[x+t]=="G")):
                    no = False
                elif(comprev[t]=="Y" and not (code[x+t]=="C" or code[x+t]=="T")):
                     no = False
                elif(comprev[t]=="M" and not (code[x+t]=="A" or code[x+t]=="C")):    
                    no = False
                elif(comprev[t]=="K" and not (code[x+t]=="G" or code[x+t]=="T")):
                    no = False
                elif(comprev[t]=="S" and not (code[x+t]=="C" or code[x+t]=="G")):
                    no = False
                elif(comprev[t]=="W" and not (code[x+t]=="A" or code[x+t]=="T")):
                    no = False
                elif(comprev[t]=="H" and not (code[x+t]=="A" or code[x+t]=="C" or code[x+t]=="T")):
                    no = False
                elif(comprev[t]=="V" and not (code[x+t]=="A" or code[x+t]=="G" or code[x+t]=="C")):
                    no = False
                elif(comprev[t]=="B" and not (code[x+t]=="C" or code[x+t]=="G" or code[x+t]=="T")):
                    no = False
                elif(comprev[t]=="D" and not (code[x+t]=="A" or code[x+t]=="G" or code[x+t]=="T")):
                    no = False
                elif(comprev[t]=="N" and not (code[x+t]=="A" or code[x+t]=="G" or code[x+t]=="T" or code[x+t]=="C")):
                    no = False


                elif (comprev[t]!=code[x+t]):
                    no = False

            if no is True:
                print ""
                print "complementary reverse"
                print z[1]
                print enzyme
                print code[:x] + start + code[x:x+lene] + end + code[x+lene:]
            else:
                no = True 







##    if(code[x : x+lene] == enzyme):
##        yes = True
##        print x
##        print("")


##  if(code[x : x+lene] == comp(enzyme)):
##      print "comp"
##      print x
      ##print code[:x] + start + code[x:x+lene] + end + code[x+lene:]
##      print("")


##  if(code[x : x+lene] == reverse(enzyme)):
##      print "reverse"
##      print x
##      print code[:x] + start + code[x:x+lene] + end + code[x+lene:]
##      print("")


##    if(code[x : x+lene] == comp(reverse(enzyme))):
##        yes = True
##        print "reverse comp"
##        print x
##        print code[:x] + start + code[x:x+lene] + end + code[x+lene:]
##        print("")


##if (yes is False):
##    print "none!"

    
