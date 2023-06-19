# Contrast mapping method
 This is a collection of notebooks for **contrast mepping** method that allows low-cost testing of field curvature and contrast of a detection objective using a Ronchi ruling, a high-frequency slide containing black and transparent lines of speficif frequency (e.g. 40 line paris/mm).

 Depending on your needs, the Ronchi ruling can be in the air (for air objectives), or immersed into an imaging medium (for immersion objectives). For our specific needs, we test **air** objectives as they image through about 20 mm of high-index medium (cleared tissues, n1.52), so we immerse the Ronchi in a chamber with oil:

![Fig2ab](https://github.com/nvladimus/lens-testing/assets/10835134/e32449eb-aa6c-4746-9a61-5c6e64b18c4b)

We are interested in field flatness, image contrast, and vignetting at the image corners, especially for CMOS camers with relatively large sensors (25 mm and higher), so that we can select the most suitable long-WD air objectives for our custom-built light-sheet microscopes.

 Objectives tested in this repository:
- PlatinumTL™ 0.9x Telecentric Lens
- Lensagon T25M-12-155I Telecentric Lens
- Mitutoyo BD Plan Apo 2x/0.055
- Mitutoyo BD Plan Apo 5x/0.14
- Mitutoyo BD Plan Apo 7.5/0.21
- Mitutoyo BD Plan Apo 10/0.28
- Mitutoyo G Plan Apo 20x/0.28(t3,5)
- Olympus MVPLAPO-1x, zoom 1x
- Olympus MVPLAPO-1x, zoom 1.25x
- Olympus MVPLAPO-1x, zoom 2x
- Olympus MVPLAPO-1x, zoom 4x
- Olympus MVPLAPO-1x, zoom 5x
- Olympus MVPLAPO-1x, zoom 6.3x
- Thorlabs Super Apochromatic 2x/0.1 TL2X-SAP
- Thorlabs Super Apochromatic 4x/0.2 TL2X-SAP
- Olympus XLFluor 4x/0.28

See [Wiki](https://github.com/nvladimus/lens-testing/wiki) for summary plots of various objectives.

The method of contrast mapping was developed by [N.Vladimirov](https://github.com/nvladimus) as a part of the [Benchtop mesoSPIM](https://github.com/mesoSPIM/benchtop-hardware) project.
