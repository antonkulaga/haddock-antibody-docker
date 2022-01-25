version development
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow haddock_antibody {
    input {
        String destination
        String name
        File antibody
        File antigen
    }

    call prepare {
        input:
            output_name = name,
            antibody = antibody,
            antigen = antigen
    }

    call files.copy as copy{
        input: files = [prepare.out], destination = destination

    }
    output {
        File results = copy.out[0]

        File ambig = prepare.ambig
        File unambig = prepare.unambig
        File antigen_active_passive = prepare.antigen_active_passive
        File antibody_active_passive = prepare.antibody_active_passive
        File antibody_pdb = prepare.antibody_pdb
        File antigen_pdb = prepare.antigen_pdb
    }

}

task prepare {

    input {
        File antibody
        File antigen
        String output_name
    }

    String antigen_name = basename(antigen, ".pdb")
    String antibody_name = basename(antibody, ".pdb")


    command {
        start.py start --antibody ~{antibody} --antigen ~{antigen} --output ~{output_name}
    }

    runtime {
        docker: "quay.io/antonkulaga/haddock-antibody@sha256:4a581a97fbf9242f23b49cab04df5a8aa3f7eaec68792159aa1e08d2e2f0541d"
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
    }

}