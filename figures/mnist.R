# Load packages -------
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load("tidyverse", "gridExtra", "png", "stringr", "magick", "colorspace", "pracma")

## MNIST Grid ----------

mnist_dir <- "../data/input/mnist/"
OUT_DIR <- "../data/output/"

dfs <- list(
  Unadjusted = read_csv(paste(c(mnist_dir, "unadjusted.csv"), collapse="")) %>% select(-Sample, -Digit),
  `Artificial noise` = read_csv(paste(c(mnist_dir, "noisy.csv"), collapse="")) %>% select(-Batch, -Sample, -Digit),
  Scale = read_csv(paste(c(mnist_dir, "noisy_scale.csv"), collapse="")) %>% select(-Batch, -Sample, -Digit),
  ComBat = read_csv(paste(c(mnist_dir, "noisy_combat.csv"), collapse="")) %>% select(-Batch, -Sample, -Digit),
  Confounded = read_csv(paste(c(mnist_dir, "noisy_confounded.csv"), collapse="")) %>% select(-Batch, -Sample, -Digit)
)


pic_rows <- c(1, 2223, 4445, 1112, 3334, 5556)

image_from_vector <- function(x) {
  return(image_read(aperm(array(as.numeric(x), c(28, 28, 1)), c(2, 1, 3))))
}

big_img <- NULL
imgs <- list()
for (name in names(dfs)) {
  img <- NULL
  df <- dfs[[name]]
  for (i in pic_rows) {
    new_img <- image_from_vector(df[i,])
    if (is.null(img)) {
      img <- new_img
    } else {
      img <- image_append(c(img, new_img))
    }
  }
  imgs[[i]] <- img
  if (is.null(big_img)) {
    big_img <- img
  } else {
    big_img <- image_append(c(big_img, img), stack = TRUE)
  }
}

# Perhaps tweak these sizes to get the labels to show up correctly.
label_width <- 65
text_background <- image_scale(image_read(array(rep(1, 28*label_width), c(28, label_width, 1))), "x100")
labels <- NULL
for (name in names(dfs)) {
  # Font size = 32
  new_label <- image_annotate(text_background, name, gravity = "center", size = 32)
  if (is.null(labels)) {
    labels <- new_label
  } else {
    # This puts the images vertically.
    labels <- image_append(c(labels, new_label), stack = TRUE)
  }
}

# I think "600" means 6x100 pixels. If you have 7 images, change it to 700.
big_img <- image_scale(big_img, "600")
labelled <- image_append(c(labels, big_img))

image_write(labelled, paste0(OUT_DIR, "/mnist.png"))
