#!/usr/bin/env python2.7

import sys
import csv
from itertools import izip
from util import *

def cvtChrom(x):
    if x.isdigit():
        return int(x)
    else:
        return x

if len(sys.argv)!=4:
    print >>sys.stderr, "Usage: joinFacetsInfo.py FACETS.SEG FACETS.SAMPLE OrigMAF"
    sys.exit()

facetsSegFile=sys.argv[1]
facetsSampFile=sys.argv[2]
origMAFFile=sys.argv[3]

facetSampDb=dict()
with open(facetsSampFile) as fp:
    cin=csv.DictReader(fp,delimiter="\t")
    for r in cin:
        facetSampDb[(r["Tumor"],r["Normal"])]=r

facetSegCol="""
ID chrom loc.start loc.end seg num.mark nhet
cnlr.median mafR segclust cnlr.median.clust mafR.clust
cf tcn lcn cf.em tcn.em lcn.em
""".strip().split()

facetSegDb=dict()
with open(facetsSegFile) as fp:
    for line in fp:
        (chrom,start,end,_,_,_,_,facetData)=line.strip().split()
        if facetData==".":
            continue
        start=int(start)+1
        key="%s:%d-%s" % (chrom,start,end)
        facetDict=dict()
        for (k,v) in zip(facetSegCol,facetData.split("|")):
            facetDict[k]=v
        (proj,tumor,normal)=facetDict["ID"].split("_")
        facetSegDb[(key,tumor,normal)]=facetDict

def computeCCFAndCopies(rr,f,purity):
    if f["lcn"]=="NA" or f["tcn"]=="NA":
        return ("NA","NA")
    purity=float(purity)
    tcn=float(f["tcn"])
    lcn=float(f["lcn"])
    ns=float(rr["t_alt_count"])
    nw=float(rr["t_ref_count"])
    Pi=1-purity
    M=tcn-lcn
    r=M
    m=lcn
    allele_fraction=ns/(ns+nw)
    if r>0 and purity>0:
        frac=min((allele_fraction*(purity*(M+m)+2*Pi)/purity)/r,1)
    else:
        #
        # If r or purity is == 0 then 
        # either no tumor cells (purity==0)
        # or no copies (r==M==tcn-lcn==0)
        return (1,0)

    rSeq=xrange(1,int(M+1))
    frac1 = [purity * r / (purity * (M+m) + 2*Pi) for r in rSeq]

    pcopies = [(f ** ns) * ((1-f)**nw) for f in frac1]

    copies=max(zip(pcopies,rSeq))[1]

    return (frac,copies)

##############################################################################
##############################################################################

events=dict()
with open(origMAFFile) as fp:
    print fp.readline(),
    cin=csv.DictReader(fp,delimiter="\t")
    for r  in cin:
        if r["Reference_Allele"]!=r["Tumor_Seq_Allele1"]:
            alt=r["Tumor_Seq_Allele1"]
        elif r["Reference_Allele"]!=r["Tumor_Seq_Allele2"]:
            alt=r["Tumor_Seq_Allele2"]
        else:
            print >>sys.stderr, "Should never get here"
            print >>sys.stderr
            print >>sys.stderr, r
            print >>sys.stderr
            sys.exit()
        pos=r["Chromosome"]+":"+r["Start_Position"]+"-"+r["End_Position"]
        tag=pos+":"+r["Reference_Allele"]+":"+alt

        tumor=r["Tumor_Sample_Barcode"].split("_")[1]
        normal=r["Matched_Norm_Sample_Barcode"].split("_")[1]
        label=tag+"::"+tumor+":"+normal

        if (tumor,normal) in facetSampDb:
            facetSampInfo=facetSampDb[(tumor,normal)]
            r["Purity"]=facetSampInfo["Purity"]
            r["Ploidy"]=facetSampInfo["Ploidy"]
            r["WGD"]="Yes" if float(facetSampInfo["dipt"]>2) else "No"

            fKey=(pos,tumor,normal)
            if fKey in facetSegDb:
                facetSegInfo=facetSegDb[(pos,tumor,normal)]
                r["seg.mean"]=facetSegInfo["cnlr.median"]
                r["cf"]=facetSegInfo["cf"]
                r["tcn"]=facetSegInfo["tcn"]
                r["lcn"]=facetSegInfo["lcn"]
                if facetSampInfo["Purity"]!="NA":
                    (est_frac,est_copies)=computeCCFAndCopies(r,facetSegInfo,float(facetSampInfo["Purity"]))
                else:
                    (est_frac,est_copies)=("NA","NA")
                r["est_frac"]=est_frac
                r["est_copies"]=est_copies
            else:
                print >>sys.stderr, fKey, "Missing in facets seg data"
                r["seg.mean"]="NA"
        else:
            print >>sys.stderr, (tumor, normal), "FACETS for this sample is missing"
            
            r["Purity"]="NULL"
            r["Ploidy"]="NULL"
            r["WGD"]="NULL"

        events[label]=r

outFields=cin.fieldnames+[  "Purity","Ploidy","WGD",
                            "seg.mean","cf","tcn","lcn",
                            "est_frac","est_copies"]
cout=csv.DictWriter(sys.stdout,outFields,delimiter="\t")
cout.writeheader()
for ki in sorted(events):
    cout.writerow(events[ki])

