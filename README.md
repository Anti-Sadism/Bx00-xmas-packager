Welcome to the 0x00-xmas-packager wiki!

# Overview

Can a person hide files inside of YouTube videos undetectably? Can YouTube serve as a secret file host? What would happen if it developed that YouTube could be used for this purpose? Steganography is almost pointless if I can only hide text. What if I could hide videos, within videos? A photo collection? Large rar files and EXE's? Anything?

If I download an innocent seeming video teeming with hidden illicit content, without knowing what's there, am I legally guilty of anything? What if it is possible to share illicit content in YouTube videos? Can we detect this?

### What kind of files can I hide?

You can hide any kind of file from images to executable code. 3GP video works best for video within video applications, .OGG files work nicely for music, and .RAR files work just dandy for photo albums, and collections of dissimilar file types. Of course, I'm just giving you suggestions on how to reduce the filesize.

### The rumored AmerisourceBergen attack in 2014:

**SkyHigh Networks** recently published information about a sophisticated attack where employees of a "certain" unnamed Fortune 500 company were sending out confidential data as YouTube videos. There wasn't anything special about these videos. Remember that the YouTube videos were not altered by the data in ways that were visually detectable?


# Academic Papers
1. (**smiling-face**) Ming Yang; Bourbakis, N., "A high bitrate information hiding algorithm for digital video content under H.264/AVC compression," in Circuits and Systems, 2005. 48th Midwest Symposium on , vol., no., pp.935-938 Vol. 2, 7-10 Aug. 2005
doi: 10.1109/MWSCAS.2005.1594256
 keywords: {multimedia communication;quantisation (signal);security of data;video codecs;video coding;1 bit;DCT coefficient block;H.264-AVC compression;advanced video codec;channel capacity;digital video content;high bitrate information hiding algorithm;lossy video codec;quantization matrix;vector quantization;video coding standard;Animation;Automatic voltage control;Bit rate;Channel capacity;Discrete cosine transforms;Image coding;Robustness;Video codecs;Video compression;Watermarking},

2. (**checkers**) Alavianmehr, M.A.; Rezaei, M.; Helfroush, M.S.; Tashk, A., "A reversible data hiding scheme for video robust against H.264/AVC compression," in Information Security and Cryptology (ISCISC), 2013 10th International ISC Conference on , vol., no., pp.1-6, 29-30 Aug. 2013
doi: 10.1109/ISCISC.2013.6767333
 keywords: {data compression;data encapsulation;image watermarking;video coding;wavelet transforms;AVC compression;H.264 compression;integer wavelet transform domain;lossless data hiding scheme;luminance component;multilevel histogram shifting mechanism;multilevel shifting mechanism;reversible data hiding scheme;secret data;uncompressed video data;video frame;watermarked image;Bit error rate;Histograms;Image coding;PSNR;Robustness;Video coding;Watermarking;H.264/AVC;data hiding;histogram shifting;lossless;multi-level;reversible;video watermarking},

#TO-DO

**Upload test videos to one of my YouTube accounts so that new users can sample the program**

Use secure keys to recover files, and have the program automatically generate one for you.

Add ffmpeg support built in to the program, so file conversions do not need to be done outside of it.

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
