cwlVersion: v1.0
class: CommandLineTool
label: PyClone
baseCommand: ["Rscript", "/home/pipeline/run_analysis_pyclone.R"]
requirements:
  - class: DockerRequirement
    dockerPull: smcheteval/pyclone:0.1

inputs:
  input_vcf:
    type: File
    inputBinding:
      position: 1

  battenberg_file:
    type: File
    inputBinding:
      position: 2

  purity_file:
    type: File
    inputBinding:
      position: 3

  num_mcmc:
    type: int
    default: 10
    inputBinding:
      position: 4

  burn_in:
    type: int
    default: 2
    inputBinding:
      position: 5

  max_snv:
    type: int
    default: 25
    inputBinding:
      position: 6

outputs:
  population:
    type: File
    outputBinding:
      glob: 1B.txt

  proportion:
    type: File
    outputBinding:
      glob: 1C.txt

  cluster_assignment:
    type: File
    outputBinding:
      glob: 2A.txt

  cocluster_assignment:
    type: File
    outputBinding:
      glob: 2B.txt
