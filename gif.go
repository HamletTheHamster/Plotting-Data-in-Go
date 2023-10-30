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
	"sort"
	"strconv"
	"sync"
)

type frameResult struct {
	Index   int
	Palette *image.Paletted
}

func main() {
	// List of PNG filenames
	pngFilenames := []string{}
	for i := 0; i < 75; i++ {
		pngFilenames = append(pngFilenames, "plots/2023-Oct-30/gifFrames/CABS"+strconv.Itoa(i)+".png")
	}

	// Use the last PNG to generate the palette
	lastPNG := openPNG(pngFilenames[len(pngFilenames)-1])
	palette := generatePalette(lastPNG)

	var frames []*image.Paletted
	var delays []int

	// Channel to receive frame results
	resultCh := make(chan frameResult, len(pngFilenames))
	var wg sync.WaitGroup

	for index, fname := range pngFilenames {
		wg.Add(1)
		go convertToPaletted(index, fname, palette, resultCh, &wg)
	}

	wg.Wait()
	close(resultCh)

	// Collect the results in a slice
	var frameResults []frameResult
	for result := range resultCh {
		frameResults = append(frameResults, result)
	}

	// Sort the frameResults by their index
	sort.Slice(frameResults, func(i, j int) bool {
		return frameResults[i].Index < frameResults[j].Index
	})

	for _, result := range frameResults {
		frames = append(frames, result.Palette)
		delays = append(delays, 5)
	}

	// Save as a GIF
	outFile, _ := os.Create("animated.gif")
	defer outFile.Close()
	gif.EncodeAll(outFile, &gif.GIF{
		Image: frames,
		Delay: delays,
	})
}

func convertToPaletted(index int, fname string, palette []color.Color, resultCh chan<- frameResult, wg *sync.WaitGroup) {
	defer wg.Done()
	img := openPNG(fname)
	palettedImage := image.NewPaletted(img.Bounds(), palette)
	draw.Draw(palettedImage, img.Bounds(), img, image.Point{}, draw.Over)
	resultCh <- frameResult{
		Index:   index,
		Palette: palettedImage,
	}
}

func generatePalette(img image.Image) []color.Color {
	paletted := image.NewPaletted(img.Bounds(), palette.Plan9)
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
