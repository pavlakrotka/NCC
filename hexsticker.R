library(hexSticker)
library(showtext)

font_add_google("Neucha", "Neucha")

s <- sticker(magick::image_read(rsvg::rsvg_raw("man/figures/trial_hexsticker.svg")),
             package="NCC", p_size=15, s_x=1, s_y=1.1, s_width=1.6, s_height=1, p_y=0.45,
             h_fill = "#F5F5F5", h_color = "#8B0000", h_size = 0.5,
             p_color = "darkred", p_family = "Neucha",
             spotlight = F,
             filename = "man/figures/NCC_hexsticker.png")

print(s)

save_sticker("man/figures/NCC_hexsticker.png", sticker = s, scale = 0.5)
