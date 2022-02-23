version development

workflow haddock_antibodies{
    input {
        String destination
        String name
        File antibody
        File antigen
        Int n_comp = 4
        Int run_number = 1
    }
}