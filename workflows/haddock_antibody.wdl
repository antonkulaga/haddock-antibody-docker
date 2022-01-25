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
    }

}

task prepare {

    input {
        File antibody
        File antigen
        String output_name
    }

    command {
        start.py start --antibody ~{antibody} --antigen ~{antigen} --output ~{output_name}
    }

    runtime {
        docker: "quay.io/antonkulaga/haddock-antibody:latest"
    }

    output {
        File out = output_name
    }

}