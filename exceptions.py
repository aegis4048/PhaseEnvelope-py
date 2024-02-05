def raise_missing_required_data_error(name, property_description):
    raise ValueError(
        "Chemical name '%s' is recognized but missing a required data (%s)." % (name, property_description))