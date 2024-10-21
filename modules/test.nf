process TEST {
    time = '1d'
    cpus = 1
    memory = '12G'

    tag "TEST"

    input:
    tuple path( peak_list )

    output:
    tuple path( "*" )

    script:
    """
    echo $peak_list >out.txt
    
    """
}