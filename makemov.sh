ffmpeg -i all%03d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p out.mp4
