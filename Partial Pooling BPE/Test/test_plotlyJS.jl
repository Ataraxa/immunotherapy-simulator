using PlotlyJS: scatter, plot, relayout!

p1 = plot(scatter(x=1:3, y=4:6))
p2 = plot(scatter(x=20:40, y=50:70))
p = [p1 p2]
relayout!(p, title_text="Side by side layout (1 x 2)")
p