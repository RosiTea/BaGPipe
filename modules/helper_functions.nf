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

    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }

}
