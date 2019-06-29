import os
configfile: 'config.yaml'

if not os.path.exists('logs/RATs'):
	os.mkdir('logs/RATs')


if config["design"]["paired"]:
	SAMPLES = list(set([ "_".join(x.split("_")[:2]) for x in FILES]))
else:
	SAMPLES = list(set([ x.rstrip(extension) for x in FILES]))


rule run_RATS:
	input:
		expand('Salmon/quants/{sample}/quant.sf', sample=SAMPLES)

	output:
		'DTU/DTU.txt'

	params:
		padj = config["DTU"]["pADJ"],
		cutoff = config["DTU"]["cutoff"],
		gtf = "Reference/reference.gtf",
		cond_control = config["design"]["condition_1"],
		cond_traitment = config["design"]["condition_2"],
		outprefix = "DTU/"

	log:
		"logs/RATs/RATs.log"

	message: ''' --- Running DTU analysis --- '''

	shell: '''
        Rscript scripts/DTU.R --cutoff={params.cutoff} \
    		--padj={params.padj} \
    		--gtf={params.gtf} \
   			--cond_controle={params.cond_control} \
   			--cond_traitment={params.cond_traitment} \
    		--outprefix={params.outprefix} \
            >{log} 2>&1 \
        '''
