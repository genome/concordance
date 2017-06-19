cwlVersion: v1.0
class: CommandLineTool
label: "concordance_tool"
baseCommand: ["python", "/opt/concordance/newConcordance.py"]
hints: 
    class: DockerRequirement
    dockerPull: "mneveau/docker-cle"
inputs: 
    reference: 
	type: File
	inputBinding: 
	    position: 3
	secondaryFiles: [".fai", "^.dict"]
    tumor_bam: 
	type: File
	inputBinding: 
	    position: 1
	secondaryFiles; [^.bai]
    normal_bam: 
	type: File
	inputBinding: 
	    position: 2
	secondaryFiles: [^.bai]
    snp: 
	type: File
	inputBinding:
	    position: 4
    output_file_name: 
	type: string
	inputBinding:  
	    position: 5
    output_genotypes_file_name: 
	type: string
	inputBinding: 
	    position: 6
outputs: 
    output_file: 
	type: File
	outputBinding:
  	    glob: "*.txt"
    output_genotypes:
	type: File
	outputBinding:
	    glob: "*.bed"


