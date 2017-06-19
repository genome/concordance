cwlVersion: v1.0
class: CommandLineTool
label: "Concordance tool"
baseCommand: ["python", "/opt/concordance/newConcordance.py"]
requirements:
    - class: DockerRequirement
      dockerPull: "mneveau/docker-cle:concordance_addition"
inputs:
    bam_1:
        type: File
        inputBinding:
            position: 1
        secondaryFiles: [.bai]
    bam_2:
        type: File
        inputBinding:
            position: 2
        secondaryFiles: [.bai]
    reference:
        type: File
        inputBinding:
            position: 3
        secondaryFiles: [.fai, "^.dict"]
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