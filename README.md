Welcome to the 0x00-xmas-packager wiki!

# Overview

### Why was xmas-Packager created?

### The Attack:

SkyHigh Networks recently published information about a sophisticated attack where employees of an unnamed Fortune 500 company were sending out confidential data as YouTube videos. There wasn't anything special about these videos. The YouTube videos were not altered by the data in ways that were visually detectable. 

In addition, the hidden data did not alter the file-size of the video. Discernibly, there was very little difference between the video and any other YouTube video, except of course, that these videos contained credit card information, source and binary code, names, passwords, and the like.

**Forensic steganography** is something that few actually study. See, it involves hiding data inside of other data, in ways that don't arouse suspicion. It's a very sophisticated science and is separate from Cryptography for its own reasons. Knowing how to use it however, can probably help you someday soon.

The attack was undetectable, that is, until the sysops noticed multiple uploads of the same video (The attackers weren't very creative even with such a sophisticated tool, and decided to use the same video to their downfall)

**Why didn't they just steal data using their own removable media?** Some corporate networks do not allow employees to use removable media. In addition, firewalls that are watched by human operators, as well as activity logs make it difficult to send a great amount (MB, GB sized) data out.

Compare this to the arrest made on an **Al-Qaeda operative** who was found to have terrorist training manuals within a pornography video in 2012. His method however, was less than sophisticated. It involved stuffing the TrueCrypt chunk within a forbidden zone of the video file, and altered the filesize of the video. Still, it is marvelous how much data was able to be fit within the file.

It goes to show you that even the most foolish of ideas can become powerful if your enemy does not consider them.

B1N@RY and his student wrote a paper titled UYR, where they demonstrate one possibility of such a powerful attack. Their attack uses LUT's and doesn't actually hide the data inside of the video. XMAS-Packager hopes to take this to the extreme where all that is needed is a single 32 bit length key and very large files can be hosted and downloaded from YouTube, even after re-encoding. It combines digital image processing with forensic steganography.

One can take advantage of this to share files without needing strong encryption provided by anonymizing services like Tor. The hiding capacity of video is fairly large. **All the feds will see, as well as your ISP, and anyone else, is legitimate viewing and sharing of YouTube videos. If you choose to just keep it on your disk, even better, you won't need to use error correction and can invisibly hide even more.** So have fun hiding your programs, pictures, and videos, folks.

So forensic steganography isn't new. Sure, there ARE programs that can do it. However, XMAS-packager is one of few video steganography programs, designed for YouTube and specifically online video websites, which can reliably use the same parameters and encode data losslessly for any video. It is not data dependent and you won't need to fool around with parameters constantly just to get something sorta-working.

XMAS-Packager uses many algorithms. But all algorithms must at least 1.) Have >30db distortion min 2.) Not alter the filesize 3.)Be antifragile, robust, or at least marginally so. If you invent any of your own, feel free to contribute.

### This Program is Dangerous to Use!

Really though, why? There's a fair amount of bit error and this is dealt with using error correcting codes. But a single bit error in a fatal program could possibly cause erratic and unpredictable behavior for a program. **PLEASE, Take care to only run programs to which the total BER (Bit Error Rate) is 0. Only run programs from trusted sources.** Don't worry though, if you follow the instructions here, you shouldn't have a single bit error to deal with.

### Watch what I mean, Let's hide a single JPEG within the first video frame.

**Original file only hidden in the first frame of a 10 second 300 frame video**
![Lena Original](http://s18.postimg.org/9y7xy1fh5/lena.jpg)

Filesize: 16KB, JPEG Lena.jpg

**File in first frame after YouTube's H.264 compression, no error code used.**
![](http://s24.postimg.org/vl6fn042d/final_result.jpg)

Filesize: 16KB, JPEG, BER = 000.29% Lena.jpg

**File in first frame hidden with error code r(50, 1) after YouTube encoding**
![](http://s16.postimg.org/6myej8xjp/lena_image.jpg)

Filesize 16KB, JPEG, BER = 0% Lena.jpg

**So using error correcting codes makes it lossless. Luckily, the same code will work for all file/video combinations.**

Follow the instructions below for using each algorithm and you'll be guaranteed to have no bit errors.

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

1. (**checkers**) Using the spatial domain PVD algorithm

2. (**smiling-face**) Using the dct 8-D vector quantization algorithm

# Academic Papers
1. (**smiling-face**) Ming Yang; Bourbakis, N., "A high bitrate information hiding algorithm for digital video content under H.264/AVC compression," in Circuits and Systems, 2005. 48th Midwest Symposium on , vol., no., pp.935-938 Vol. 2, 7-10 Aug. 2005
doi: 10.1109/MWSCAS.2005.1594256
 keywords: {multimedia communication;quantisation (signal);security of data;video codecs;video coding;1 bit;DCT coefficient block;H.264-AVC compression;advanced video codec;channel capacity;digital video content;high bitrate information hiding algorithm;lossy video codec;quantization matrix;vector quantization;video coding standard;Animation;Automatic voltage control;Bit rate;Channel capacity;Discrete cosine transforms;Image coding;Robustness;Video codecs;Video compression;Watermarking},

2. (**checkers**) Alavianmehr, M.A.; Rezaei, M.; Helfroush, M.S.; Tashk, A., "A reversible data hiding scheme for video robust against H.264/AVC compression," in Information Security and Cryptology (ISCISC), 2013 10th International ISC Conference on , vol., no., pp.1-6, 29-30 Aug. 2013
doi: 10.1109/ISCISC.2013.6767333
 keywords: {data compression;data encapsulation;image watermarking;video coding;wavelet transforms;AVC compression;H.264 compression;integer wavelet transform domain;lossless data hiding scheme;luminance component;multilevel histogram shifting mechanism;multilevel shifting mechanism;reversible data hiding scheme;secret data;uncompressed video data;video frame;watermarked image;Bit error rate;Histograms;Image coding;PSNR;Robustness;Video coding;Watermarking;H.264/AVC;data hiding;histogram shifting;lossless;multi-level;reversible;video watermarking},
