Welcome to the 0x00-xmas-packager wiki!

# Overview

Can a person hide files inside of YouTube videos undetectably? Can YouTube serve as a secret file host? What would happen if it developed that YouTube could be used for this purpose? Steganography is almost pointless if I can only hide text. What if I could hide videos, within videos? A photo collection? Large rar files and EXE's? Anything?

If I download an innocent seeming video teeming with hidden illicit content, without knowing what's there, am I legally guilty of anything? What if it is possible to share illicit content in YouTube videos? Can we detect this?

### What kind of files can I hide?

You can hide any kind of file from images to executable code. 3GP video works best for video within video applications, .OGG files work nicely for music, and .RAR files work just dandy for photo albums, and collections of dissimilar file types. Of course, I'm just giving you suggestions on how to reduce the filesize.

### The rumored AmerisourceBergen attack in 2014:

**SkyHigh Networks** recently published information about a sophisticated attack where employees of a "certain" unnamed Fortune 500 company were sending out confidential data as YouTube videos. There wasn't anything special about these videos. Remember that the YouTube videos were not altered by the data in ways that were visually detectable?


**Original file only hidden in the first frame of a 10 second 300 frame video**
![Lena Original](http://s18.postimg.org/9y7xy1fh5/lena.jpg)

Filesize: 16KB, JPEG Lena.jpg

**Recovered File in first frame after YouTube's H.264, no error code used.**
![](http://s24.postimg.org/vl6fn042d/final_result.jpg)

Filesize: 16KB, JPEG, BER = 000.29% Lena.jpg

**Recovered File in first frame hidden with error code r(50, 1) after YouTube**
![](http://s16.postimg.org/6myej8xjp/lena_image.jpg)

# Some Forensic Steganography BG:

1. _**Robustness** - Indifferent to attacks. Hidden content remains after corruption of cover media._

2. _**Fragility** - Breaks when attacked. Hidden content degrades after corruption of cover media._

3. _**Antifragility** - Strengthens with attacks (To a certain threshold). Hidden content becomes better defined as it is attacked._

4. _**Attack**- Anything which is performed on the stegano file. JPEG compression, YouTube compression, scaling, skewing, and the like._

# Anti-Fragility

Many of the algorithms used in xmas-packager are fragile or robust. This is a problem. We hope to eventually devise an antifragile algorithm. 

Following statistician **Nassim Nicholas Taleb's** wide-selling work, _Anti-Fragile_, we hope to eventually devise a way to make _anti-fragile_ embedding attempts within video. Taleb mentions that many areas of science can stand to benefit from a careful formulation of Fragility. I believe forensic steganography is one of those areas. Witnessed above is an example of the well known "Avalanche Effect" in computer science, where a single bit error can potentially trigger more bit errors within a file and the entire future run length of the file is corrupted from a single bit; this file is fragile.

Anti-Fragility is strongly desired for steganography and the person who finds an algorithm will surely be hailed as the champion of steganography (Possibly with his or her face on a box of Wheaties :) ). 

The antifragile algorithm is truly the breakfast of champions. For example, PVD algorithms are antifragile to contrast filters, and contrast filters are lossy information-wise. A video steganography algorithm anti-fragile to YouTube compression would cause data to be better defined each time the video is encoded.

### Current Strategy

Error correcting codes, such as the **r(n, 1)** repetition code help, (where n is the length of the codeword) but error correcting codes themselves can be quite fragile. That is to say, as the number of error bits increase, the closer the message comes to being erroneously corrected as a false positive.

Imagine I have a r (8, 1) codeword with the message "0" and codeword "00000000"; after sending it to YouTube I get "01011011". There are more 1's than zeros and if this is the r(8, 1) code, then I'll correct it, incorrectly, as a 1 bit. This is the fragility of error correcting codes.

# Performance

Certain algorithms perform better than others. If you are hoping to send data to YouTube, then use (checkers). Though (smiling-face) can be used, it is low-bitrate. We recommend it for disk storage. But please read (here) on how to balance invisibility with robustness. Contributions are very much welcome. 

Needed work involves producing more algorithms and testing known ones against **forensic analysis**. It is known that (checkers) can be detected with a chi-squared attack. (smiling-face) is more difficult to detect however.

# Getting Started 

With a minimalist style:

QUICKSTART: https://github.com/Anti-Government/0x00-xmas-packager/wiki/Quick-Start

1. (**checkers**) Using the spatial domain PVD algorithm [here](https://github.com/Anti-Government/0x00-xmas-packager/wiki/Algorithm:-(checkers))

2. (**smiling-face**) Using the Discrete Cosine Transform 8-D vector quantization algorithm [here](https://github.com/Anti-Government/0x00-xmas-packager/wiki/Algorithm:-(smiling-face))

# Academic Papers
1. (**smiling-face**) Ming Yang; Bourbakis, N., "A high bitrate information hiding algorithm for digital video content under H.264/AVC compression," in Circuits and Systems, 2005. 48th Midwest Symposium on , vol., no., pp.935-938 Vol. 2, 7-10 Aug. 2005
doi: 10.1109/MWSCAS.2005.1594256
 keywords: {multimedia communication;quantisation (signal);security of data;video codecs;video coding;1 bit;DCT coefficient block;H.264-AVC compression;advanced video codec;channel capacity;digital video content;high bitrate information hiding algorithm;lossy video codec;quantization matrix;vector quantization;video coding standard;Animation;Automatic voltage control;Bit rate;Channel capacity;Discrete cosine transforms;Image coding;Robustness;Video codecs;Video compression;Watermarking},

2. (**checkers**) Alavianmehr, M.A.; Rezaei, M.; Helfroush, M.S.; Tashk, A., "A reversible data hiding scheme for video robust against H.264/AVC compression," in Information Security and Cryptology (ISCISC), 2013 10th International ISC Conference on , vol., no., pp.1-6, 29-30 Aug. 2013
doi: 10.1109/ISCISC.2013.6767333
 keywords: {data compression;data encapsulation;image watermarking;video coding;wavelet transforms;AVC compression;H.264 compression;integer wavelet transform domain;lossless data hiding scheme;luminance component;multilevel histogram shifting mechanism;multilevel shifting mechanism;reversible data hiding scheme;secret data;uncompressed video data;video frame;watermarked image;Bit error rate;Histograms;Image coding;PSNR;Robustness;Video coding;Watermarking;H.264/AVC;data hiding;histogram shifting;lossless;multi-level;reversible;video watermarking},

#TO-DO

Add feature for (checkers) that allows for assume pepper, and feature to randomly decide salt or pepper when pixel difference is beyond 2T+G.

**Upload test videos to one of my YouTube accounts so that new users can sample the program**

Use secure keys to recover files, and have the program automatically generate one for you.

Add support for strong encryption taking into account steganalysis detection.

Add ffmpeg support built in to the program, so file conversions do not need to be done.

Speed up the 2D DCT with a faster algorithm.

Add support for calculating the PSNR and MSE.

**Use PSNR and BER to check if the data is reliably hidden before uploading to YouTube.**

**Create a better detection utlity to detect steganography in videos (Of course, I want to break my own program!)**

**Ability to specify the YouTube URL directly, without having to manually download anything and convert with ffmpeg**

Integrate it with live webcam feed. To do so, devise a computationally efficient algorithm. The DCT is slow.

# Other Ideas

Use the Discrete Chirplet Transform for the video embedding. It has shown to have supreme robustness over the DCT and is antifragile to blur and rotate operations.

Add a seeker which can seek a video that has data hidden in it. Allow the user only to embed within certain sections of the video, rather than the whole slice.

Add audio steganography 

Allow combinations of algorithms to be performed at once. Some work well with others.

Allow the setting of BCH Syndrome Code and Hamming/CRC code instead of Repetition Error correcting code which would allow for greater payload per codeword length
