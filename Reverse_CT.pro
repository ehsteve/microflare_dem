PRO Reverse_CT

; This program reverses the colors in the current color table.

TVLCT, r, g, b, /Get
TVLCT, Reverse(r), Reverse(g), Reverse(b)
END