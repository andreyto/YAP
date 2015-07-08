import argh
from argh import arg
import os, shutil, glob
from subprocess import check_call


@arg("mothur_path", type=os.path.abspath, help="Path to Mothur binary")
@arg("pcr_target_fasta",type=os.path.abspath,
        help="Path to multi-fasta that contains PCR target reference sequence (e.g. ecoli)")
@arg("pcr_target_idseq", 
        help="Sequence ID of PCR target reference sequence in pcr_target_fasta file")
@arg("--pcr-target-padding",
        help="When cutting PCR target reference region by primers, expand by this many nucleotides on both sides")
@arg("primer_forward", 
        help="String with the forward primer sequence")
@arg("primer_reverse", 
        help="String with the reverse primer sequence")

def call_mothur_cmds(mothur_path,cmds,*l,**kw):
    check_call([mothur_path,'#'+';'.join(cmds)],*l,**kw)

def find_mothur_output(patt):
    found = list(glob.glob(patt))
    assert len(found) == 1
    return found[0]

def extract_template_ali_range(template_name,ali_report_file):
    import csv
    start = 10**10
    end = -1
    template_length = 0
    template_found = False
    with open(ali_report_file,"r") as inp_file:
        inp = csv.DictReader(inp_file,delimiter="\t")
        for rec in inp:
            if template_name is None:
                template_name = rec["TemplateName"]
            if template_name == rec["TemplateName"]:
                start = min(start,int(rec["TemplateStart"]))
                end = max(end,int(rec["TemplateEnd"]))
                template_length = rec["TemplateLength"]
                template_found = True

    assert template_found
    return dict(start=start,end=end,template_length=template_length)

def pcr_reference_alignment(mothur_path,
        pcr_target_fasta,
        pcr_target_idseq,
        primer_forward,
        primer_reverse,
        reference_ali,
        pcr_target_padding=50,
        ):
    
    scratch="scratch"
    if os.path.exists(scratch):
        shutil.rmtree(scratch)
    os.makedirs(scratch)
    
    with open(os.path.join(scratch,"pcr_target.accnos"),"w") as out:
        out.write("{}\n".format(pcr_target_idseq))

    with open(os.path.join(scratch,"input.oligos"),"w") as out:
            out.write("primer {} {}\n".format(
                primer_forward,
                primer_reverse)
                )
    
    shutil.copy(pcr_target_fasta,os.path.join(scratch,"pcr_target.fasta"))
    shutil.copy(reference_ali,os.path.join(scratch,"ref.align"))
    call_mothur_cmds(mothur_path,[
        "get.seqs(fasta=pcr_target.fasta,accnos=pcr_target.accnos)",
        "pcr.seqs(fasta=current,oligos=input.oligos,keepprimer=F)",
        "align.seqs(fasta=current,reference=pcr_target.pick.fasta)"
        ],
        cwd=scratch)
    ali_report_file = find_mothur_output(os.path.join(scratch,"*.align.report"))
    tpl_range = extract_template_ali_range(template_name=pcr_target_idseq,ali_report_file=ali_report_file)
    start = max(1,tpl_range["start"]-pcr_target_padding)
    end = min(tpl_range["template_length"],tpl_range["end"]+pcr_target_padding)
    call_mothur_cmds(mothur_path,[
        "pcr.seqs(fasta=pcr_target.pick.fasta,start={},end={})".format(start,end),
        "align.seqs(fasta=current,reference=ref.align)",
        "pcr.seqs(fasta=ref.align,ecoli=pcr_target.pick.pcr.align,keepdots=F)"
        ],
        cwd=scratch)
    shutil.move(os.path.join(scratch,"ref.pcr.align"),".")


if __name__ == "__main__":
    argh.dispatch_commands([pcr_reference_alignment])

