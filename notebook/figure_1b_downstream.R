# Redraws manuscript Figure 1b with downstream DOC as the outcome of the
# structural causal model (model b1_downstream, notebook/scm_rebuild.qmd).
# Solid skyblue = positive effect with P > 0.90; solid coral = negative with
# P > 0.90; dashed grey = no-evidence (P <= 0.90 in both directions).
# Numbers are posterior mean standardized effects; the upstream-DOC path is
# the inheritance slope in mg per mg.
suppressMessages({library(ggplot2); library(here)})

# ---- node positions (x, y on a 0-10 canvas) --------------------------------
nodes <- data.frame(
  name = c("Max dam\nheight", "Number of\ndams", "Channel\ngradient", "Water\nresidence time",
           "Solar\nradiation", "Macrophytes", "Phytoplankton", "Soil litter\ncover",
           "Upstream\nDOC", "Downstream\nDOC"),
  x = c(0.9, 0.9, 2.9, 3.2, 4.6, 6.6, 6.2, 6.2, 8.8, 8.8),
  y = c(6.6, 4.6, 2.0, 5.0, 9.3, 7.2, 5.0, 2.9, 8.8, 5.4))
rownames(nodes) <- c("dh", "nd", "gr", "wrt", "sol", "mac", "pla", "lit", "up", "out")

# ---- edges: from, to, label, class, curvature, label offset ----------------
edges <- read.csv(text = 'from,to,label,class,curv,lx,ly
dh,wrt,,none,0,0,0
nd,wrt,0.06,pos,0,-0.1,0.35
gr,wrt,-0.08,neg,0,0.45,0
wrt,pla,0.11,pos,0,0,0.30
wrt,mac,,none,0,0,0
wrt,out,,none,0.18,0,0
sol,mac,0.14,pos,-0.15,0.55,0.55
sol,pla,0.19,pos,0.30,-0.85,0.10
sol,lit,,none,0.30,0,0
mac,out,0.08,pos,0,0,0.32
pla,out,0.14,pos,0,0,0.30
lit,out,0.07,pos,0,-0.2,0.32
up,out,0.95,pos,0,0.45,0')

# shorten each arrow so it starts/ends at the box border, not the centre
a <- merge(edges, nodes, by.x = "from", by.y = "row.names")
a <- merge(a, nodes, by.x = "to", by.y = "row.names", suffixes = c("", "2"))
dx <- a$x2 - a$x; dy <- a$y2 - a$y; L <- sqrt(dx^2 + dy^2)
a$x0 <- a$x + dx / L * 0.85; a$y0 <- a$y + dy / L * 0.55
a$x1 <- a$x2 - dx / L * 0.90; a$y1 <- a$y2 - dy / L * 0.55
a$xm <- (a$x0 + a$x1) / 2 + a$lx; a$ym <- (a$y0 + a$y1) / 2 + a$ly

pal <- c(pos = "#3B9AD1", neg = "#E8603C", none = "grey55")   # skyblue / coral / grey
lty <- c(pos = "solid",   neg = "solid",   none = "dashed")
arr <- arrow(length = unit(0.16, "cm"), type = "closed")

p <- ggplot()
# straight and curved arrows drawn separately (geom_curve takes one curvature at a time)
st <- subset(a, curv == 0)
p <- p + geom_segment(data = st, aes(x = x0, y = y0, xend = x1, yend = y1),
                      colour = pal[st$class], linetype = lty[st$class], linewidth = 0.7, arrow = arr)
for (i in which(a$curv != 0))
  p <- p + geom_curve(data = a[i, ], aes(x = x0, y = y0, xend = x1, yend = y1),
                      curvature = a$curv[i], colour = pal[a$class[i]],
                      linetype = lty[a$class[i]], linewidth = 0.7, arrow = arr)
p <- p +
  geom_text(data = subset(a, label != ""), aes(x = xm, y = ym, label = label), size = 4.2) +
  geom_label(data = nodes[rownames(nodes) != "out", ], aes(x = x, y = y, label = name),
             size = 4.1, label.padding = unit(0.3, "lines"), label.r = unit(0.05, "lines")) +
  geom_label(data = nodes["out", ], aes(x = x, y = y, label = name), size = 4.4,
             fontface = "bold", label.padding = unit(0.45, "lines"),
             label.r = unit(0.9, "lines"), linewidth = 0.9, colour = "#C0208E") +
  annotate("text", x = 0.15, y = 9.8, label = "(b)", size = 6.5, hjust = 0) +
  xlim(0, 10) + ylim(1.2, 10) + theme_void()
ggsave(here("images", "Figure_1b_downstream.png"), p, width = 9, height = 5.6, dpi = 300, bg = "white")
