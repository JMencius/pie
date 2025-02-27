import logging


def check_filters(filters: dict) -> bool:
    #filters = {"only_snv": only_snv, "only_indel": only_indel, "only_sv": only_sv, "no_snv": no_snv, "no_indel": no_indel, "no_sv": no_sv, "no_double": no_double}
    only_count = ([filters["only_snv"], filters["only_indel"], filters["only_sv"]]).count(True)
    if only_count > 1:
        logging.error("Only one of --only_snv, --only_indel, or --only_sv can be set at a time.")
        return False
        
    if filters["only_snv"] and filters["no_snv"]:
        logging.error("Can not set --only_snv and --no_snv simultaneously")
        return False

    if filters["only_indel"] and filters["no_indel"]:
        logging.error("Can not set --only_indel and --no_indel simultaneously")
        return False

    if filters["only_sv"] and filters["no_sv"]:
        logging.error("Can not set --only_sv and --no_sv simultaneously")
        return False
    
    return True
