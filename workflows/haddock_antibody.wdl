version development
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow haddock_antibody {
    input {
        String destination
        String name
        File antibody
        File antigen
        Int n_comp = 4
        Int run_number = 1
        String scheme = "c"
        Boolean heavy_only = false
        Float antibody_cutoff = 0.0
    }

    call prepare {
        input:
            output_name = name,
            antibody = antibody,
            antigen = antigen,
            scheme = scheme,
            destination = destination,
            heavy_only = heavy_only,
            antibody_cutoff = antibody_cutoff
    }

    call files.copy as copy_antibody{
        input: files = [
        prepare.antibody_pdb,
        prepare.antigen_pdb,
        prepare.ambig,
        prepare.unambig,
        prepare.run_param,
        prepare.antigen_active_passive,
        prepare.antibody_active_passive,
    ], destination = destination + "/" + name
    }

    output {
        File results = copy_antibody.destination_folder
        File antibody_pdb =  copy_antibody.out[0]
        File antigen_pdb =  copy_antibody.out[1]
        File ambig = copy_antibody.out[2]
        File unambig = copy_antibody.out[3]
        File run_param = copy_antibody.out[4]
        File antigen_active_passive = copy_antibody.out[5]
        File antibody_active_passive = copy_antibody.out[6]
    }

}

task prepare {

    input {
        File antibody
        File antigen
        String output_name
        String destination
        String scheme = "c"
        Boolean heavy_only
        Float antibody_cutoff = 0.0
    }

    String antigen_name = basename(antigen, ".pdb")
    String antibody_name = basename(antibody, ".pdb")


    command {
        start.py start --scheme ~{scheme} --cutoff ~{antibody_cutoff} --antibody ~{antibody} --antigen ~{antigen} ~{if(heavy_only) then "--heavy_only true" else ""} --output ~{output_name} --project ~{destination}/~{output_name}
    }

    runtime {
        docker: "quay.io/antonkulaga/haddock-antibody"
        shell: "/usr/local/bin/_entrypoint.sh"
    }

    output {
        File out = output_name
        File ambig = output_name + "/" + "antibody-antigen-ambig.tbl"
        File unambig = output_name + "/" + "antibody-unambig.tbl"
        File antigen_active_passive = output_name + "/" + antigen_name + "_tidy_antigen_active_passive.txt"
        File antibody_active_passive = output_name + "/" + antibody_name + "_antibody_active_passive.txt"
        File antibody_pdb = output_name + "/" + antibody_name + "_HADDOCK_tidy.pdb"
        File antigen_pdb = output_name + "/" + antigen_name + "_tidy.pdb"
        File run_param = output_name + "/" + "run.param"
    }

}