# adjuster_plot_utils.R
# Common utility functions for adjuster plotting

format_adjuster_label <- function(adjuster_name) {
  if (adjuster_name == "unadjusted") return("Unadjusted")
  if (adjuster_name == "combat") return("ComBat")
  if (adjuster_name == "combat_sup") return("ComBat-Sup")
  if (adjuster_name == "combat_mean") return("Combat_mean")
  if (adjuster_name == "mnn") return("MNN")
  if (adjuster_name == "mnn_centered") return("Mnn_centered")
  return(tools::toTitleCase(adjuster_name))
}

get_classifier_ordering <- function(mxe_data, classifier_name) {
  mxe_data %>%
    filter(classifier_label == classifier_name & !is.na(classifier_label)) %>%
    group_by(adjuster) %>%
    summarise(mean_mcc = mean(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_mcc))
}

create_adjuster_labels <- function(adjuster_order) {
  sapply(adjuster_order$adjuster, format_adjuster_label)
}
