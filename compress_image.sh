#!/bin/bash

convert TennessenVA.tif -set colorspace RGB -layers flatten -alpha off -compress lzw -depth 8 \
-density 600 -adaptive-resize 4500x2400  TennessenVA_compressed.tif

convert TennessenVA_compressed.tif TennessenVA_compressed.pdf

convert TennessenStats.tif -set colorspace RGB -layers flatten -alpha off -compress lzw -depth 8 \
-density 600 -adaptive-resize 4500x2400  TennessenStats_compressed.tif

convert TennessenStats_compressed.tif TennessenStats_compressed.pdf


convert AllLoads.tif -set colorspace RGB -layers flatten -alpha off -compress lzw -depth 8 \
-density 600 -adaptive-resize 4500x2400  AllLoads_compressed.tif

convert AllLoads_compressed.tif AllLoads_compressed.pdf

convert BurdenRatioLoad.tif -set colorspace RGB -layers flatten -alpha off -compress lzw -depth 8 \
-density 600 -adaptive-resize 4500x2400  BurdenRatioLoad_compressed.tif

convert BurdenRatioLoad_compressed.tif BurdenRatioLoad_compressed.pdf

convert AllMutations.tif -set colorspace RGB -layers flatten -alpha off -compress lzw -depth 8 \
-density 600 -adaptive-resize 4500x2400  AllMutations_compressed.tif

convert AllMutations_compressed.tif AllMutations_compressed.pdf

convert BurdenRatioMutations.tif -set colorspace RGB -layers flatten -alpha off -compress lzw -depth 8 \
-density 600 -adaptive-resize 4500x2400  BurdenRatioMutations_compressed.tif

convert BurdenRatioMutations_compressed.tif BurdenRatioMutations_compressed.pdf