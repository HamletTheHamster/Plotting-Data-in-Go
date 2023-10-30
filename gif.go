package main

import (
	"image"
  "image/color"
	"image/color/palette"
	"image/draw"
	"image/gif"
	"image/png"
	"log"
	"os"
  "fmt"
  "strconv"
)

func main() {
	// List of PNG filenames
	pngFilenames := []string{}
	for i := 0; i < 5; i++ {
		pngFilenames = append(pngFilenames, "plots/2023-Oct-30/gifFrames/CABS"+strconv.Itoa(i)+".png")
	}

	// Use the last PNG to generate the palette
	lastPNG := openPNG(pngFilenames[len(pngFilenames)-1])
	palette := generatePalette(lastPNG)

	var frames []*image.Paletted
	var delays []int

	for _, fname := range pngFilenames {
    fmt.Printf(".")
		img := openPNG(fname)
		palettedImage := image.NewPaletted(img.Bounds(), palette)
		draw.Draw(palettedImage, img.Bounds(), img, image.Point{}, draw.Over)
		frames = append(frames, palettedImage)
		delays = append(delays, 50)
	}

	// Save as a GIF
	outFile, _ := os.Create("animated.gif")
	defer outFile.Close()
	gif.EncodeAll(outFile, &gif.GIF{
		Image: frames,
		Delay: delays,
	})
}

func generatePalette(img image.Image) []color.Color {
	paletted := image.NewPaletted(img.Bounds(), palette.Plan9) // Use Plan9 as a base
	draw.Draw(paletted, img.Bounds(), img, image.Point{}, draw.Over)
	return paletted.Palette
}

func openPNG(fname string) image.Image {
	f, err := os.Open(fname)
	if err != nil {
		log.Fatalf("Failed to open image: %s", err)
	}
	defer f.Close()

	img, err := png.Decode(f)
	if err != nil {
		log.Fatalf("Failed to decode image: %s", err)
	}
	return img
}
