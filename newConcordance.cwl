cwlVersion: v1.0
class: CommandLineTool
label: "concordance_tool"
baseCommand: 
#which Docker container?
hints:
    class: DockerRequirement
    dockerPull:
arguments:
inputs:
    reference:
	type: File
	inputBinding:
	    prefix: "-R"
	    position: 1
	secondaryFiles: [".fai", "^.dict"]
    tumor_bam:
	type: File
	inputBinding:
	    prefix: "-I:tumor"
	    position: 2
	secondaryFiles; [^.bai]
    normal_bam:
	type: File
	inputBinding:
	    preifix: "-I:normal"
	    position: 3
	secondaryFiles: [^.bai]
    snp:
	type: File
	inputBinding:
	    prefix: "-S"
	    position: 4
outputs:
    output_file:
	type: File
	outputBinding:
  	    glob: "output.bed"
    output_genotypes:
	type: File
	outputBinding:
	    glob: "genotype_output.txt"


