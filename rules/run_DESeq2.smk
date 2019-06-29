import os
configfile: "config.yaml"

if not os.path.exists('logs/DESeq2'):
	os.mkdir('logs/DESeq2')

rule DESeq2:
	input:
		xpdesign = "experimentalDesign.txt",
		RPKM = "DEG/RPKM.txt"

	output:
		"DEG/genes_up.txt",
		"DEG/genes_down.txt",
		"DEG/tair_ids.txt"

	params:
		padj = config["DEG"]["pADJ"],
		lfc = config["DEG"]["LFC"],
		description = "Reference/description.txt",
		outprefix = "DEG/"

	log:
		"logs/DESeq2/DESeq2.log"

	message: ''' --- DEG analysis (DESeq2) --- '''

	shell: '''
        Rscript scripts/DEG.R --padj={params.padj} \
            --lfc={params.lfc} \
            --xpdesign={input.xpdesign} \
            --description={params.description} \
            --outprefix={params.outprefix} \
            >{log} 2>&1
        '''
