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
    }

    call prepare {
        input:
            output_name = name,
            antibody = antibody,
            antigen = antigen
    }

    call files.copy as copy_antibody{
        input: files = [
        prepare.antibody_pdb,
        prepare.antigen_pdb,
        prepare.ambig,
        prepare.unambig,
        prepare.antigen_active_passive,
        prepare.antibody_active_passive
    ], destination = destination + "/" + name
    }

    call run_param as params {
        input:
            antibody = copy_antibody.out[0],
            antigen = copy_antibody.out[1],
            ambig = copy_antibody.out[2],
            unambig  = copy_antibody.out[3],
            project = copy_antibody.destination_folder,
            n_comp = n_comp,
            run_number = run_number,
            output_name = name
    }

    call files.copy as copy_param {
        input: files = [params.out],
            destination = destination + "/" + name
            }

    output {
        File results = copy_antibody.destination_folder
        File antibody_pdb =  copy_antibody.out[0]
        File antigen_pdb =  copy_antibody.out[1]
        File ambig = copy_antibody.out[2]
        File unambig = copy_antibody.out[3]
        File antigen_active_passive = prepare.antigen_active_passive
        File antibody_active_passive = prepare.antibody_active_passive
        File run_param = copy_param.out[0]
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

task run_param {
    input {
        File antibody
        File antigen
        File ambig
        File unambig
        String project
        String? haddock_dir
        Int n_comp = 2
        Int run_number = 1
        String output_name
    }

    command {
        start.py run_param --antibody ~{antibody} --antigen ~{antigen} \
        --ambig ~{ambig} --unambig ~{unambig} --project ~{project} \
        ~{"--haddock_dir" +haddock_dir} --n_comp ~{n_comp} --run_number ~{run_number}
    }

    runtime {
        docker: "quay.io/antonkulaga/haddock-antibody@sha256:4a581a97fbf9242f23b49cab04df5a8aa3f7eaec68792159aa1e08d2e2f0541d"
        shell: "/usr/local/bin/_entrypoint.sh"
    }

    output {
        File out = output_name + "/" + "run.param"
    }
}