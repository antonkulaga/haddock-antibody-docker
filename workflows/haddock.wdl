version development

workflow haddock {
    input {
        File project
        File antibody
        Int run_number = 1
        Int structures = 10000
        Int refined_structures = 400
        Int analysed_structures = 400
        Boolean verbose = true
    }

    String project_path = project


    call prepare_run{
        input:
            project_path = project_path,
            antibody = antibody,
            run_number = run_number,
    }

    call patch_run {
        input:
            run = prepare_run.run,
            project_path = project_path,
            run_number = run_number,
            antibody = antibody,
            structures = structures,
            refined_structures = refined_structures,
            analysed_structures = analysed_structures,
            verbose = verbose
   }

    call run {
        input:
            project_path = project_path,
            run = patch_run.out
    }

    output {
        File out = project
        File run_cns = patch_run.run_cns
        File run_folder = run.out
    }


}

task prepare_run{
    input {
        String project_path
        File antibody
        Int run_number = 1
    }

    String run_name = "run"+run_number

    command {
        cd ~{project_path}
        python $HADDOCK/Haddock/RunHaddock.py
    }

    runtime {
        docker: "quay.io/comp-bio-aging/haddock2:latest"
        shell: "/usr/local/bin/_entrypoint.sh"
        docker_volume1: project_path + ":" + project_path
    }

    output {
        File run = project_path + "/" + run_name
    }
}

task patch_run {
    input {
        File run
        String project_path
        Int run_number
        File antibody
        Int structures
        Int refined_structures
        Int analysed_structures
        Boolean verbose = true
    }
    String run_name = "run"+run_number

    command {
        cd ~{project_path}
        molprobity.py ~{basename(antibody)}
        echo "checked molprobity"
        cd ~{run_name}
        patch.py run.cns --structures ~{structures} --refined ~{refined_structures} --analysed ~{analysed_structures} --verbose ~{verbose}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/haddock2:tools"
        shell: "/usr/local/bin/_entrypoint.sh"
        docker_volume1: project_path + ":" + project_path
    }

    output  {
        File out = project_path + "/" + run_name
        File run_cns = project_path + "/" + run_name + "/" + "run.cns"
    }
}

task run {
    input {
        String project_path
        File run #mostly for dependency reasons
    }

    command {
        cd ~{project_path}/~{basename(run)}
        python $HADDOCK/Haddock/RunHaddock.py
    }

    runtime {
        docker: "quay.io/comp-bio-aging/haddock2:latest"
        shell: "/usr/local/bin/_entrypoint.sh"
        docker_volume1: project_path + ":" + project_path
    }

    output  {
        File out = run
    }
}