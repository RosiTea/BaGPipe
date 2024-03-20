def validate_parameters() {
    def errors = 0

    if (params.manifest) {
        manifest=file(params.manifest)
        if (!manifest.exists()) {
            log.error("The manifest file specified does not exist.")
            errors += 1
        }
    } else {
        log.error("No manifest file specified. Please specify one using the --manifest option.")
        errors += 1
    }

    if (params.genotype_method != "unitig" && params.genotype_method != "pa" && params.genotype_method != "snp") {
	log.error("Invalid genotype method. Please use one of the three options: unitig|pa|snp .")
	errors += 1
    }

    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }

}
