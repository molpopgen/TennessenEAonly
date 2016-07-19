#!/bin/bash

convert TennessenVA.tif -set colorspace RGB -layers flatten -alpha off -compress lzw -depth 8 \
-density 600 -adaptive-resize 4500x2400  TennessenVA_compressed.tif

convert TennessenVA_compressed.tif TennessenVA_compressed.pdf
