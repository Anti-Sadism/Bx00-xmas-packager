Welcome to the 0x00-xmas-packager wiki!

# Overview

### What is xmas-Packager?

Hide any file(s) inside of YouTube videos. Create "two-sided" videos, where "one-side" has hidden videos and media or file content.

### Why was xmas-Packager created?

Can a person hide files inside of YouTube videos undetectably? Can YouTube serve as a secret file host? What would happen if it developed that YouTube could be used for this purpose? Steganography is almost pointless if I can only hide text. What if I could hide videos, within videos? A photo collection? Large rar files and EXE's? Anything?

If I download an innocent seeming video teeming with hidden illicit content, without knowing what's there, am I legally guilty of anything? What if it is possible to share illicit content in YouTube videos? Can we detect this?

### Are the files hidden in the YouTube videos discoverable?

I really want to say yes, actually.

YES, but not cleanly. They need to guess the parameters. General attacks, such as the chi-squared attack may work, but YouTube's compression reduces detection because the forenic analyzer needs to guess which error correcting code you're using.

We're waiting for someone to come up with a "clean" detection strategy, but visual attacks are difficult for some algorithms, where for others, steganography is dead simple to spot, and the name of the algorithm is a joke about this. (checkers) A current area of focus is to devise more invisible algorithms which show high robustness. One future idea is to use the Discrete Chirplet transform- it shows promise for this; with its superior robustness that is nearly antifragile to some image compression operations.

### What kind of files can I hide?

You can hide any kind of file from images to executable code. 3GP video works best for video within video applications, .OGG files work nicely for music, and .RAR files work just dandy for photo albums, and collections of dissimilar file types. Of course, I'm just giving you suggestions on how to reduce the filesize.

### The rumored AmerisourceBergen attack in 2014:

**SkyHigh Networks** recently published information about a sophisticated attack where employees of a "certain" unnamed Fortune 500 company were sending out confidential data as YouTube videos. There wasn't anything special about these videos. Remember that the YouTube videos were not altered by the data in ways that were visually detectable?

![Attack](http://www.tripwire.com/state-of-security/wp-content/uploads/2014/11/videodata.png)

In addition, the hidden data did not alter the file-size of the video. Discernibly, there was very little difference between the video and any other YouTube video, except of course, that these videos may have contained **credit card information**, **health records** **source and binary computer code**, **names**, **ID**, **passwords**, and the like. 

A personal source "somehow" familiar with the incident ;) emailed us that the attackers used **OpenPuff** but... it's known that OpenPuff's video steganography isn't robust enough against YouTube's encoding. Actually, there don't currently seem to be any programs reliably robust to it. Even those which attempt to be have bit error and the files aren't the same after encoding. **Guess the attackers must've written their own binaries?**

**Forensic steganography** is something that few actually study. See, it involves hiding data inside of other data, in ways that don't arouse suspicion. It's a very sophisticated science and is separate from Cryptography for its own reasons. Knowing how to use it however, can probably help you someday soon.

The attack was undetectable, that is, until the sysops noticed multiple uploads of the same video (The attackers weren't very creative even with such a sophisticated tool, and decided to use the same video to their downfall)

**Why didn't they just steal data using their own removable media?** Some corporate networks do not allow employees to use removable media. In addition, firewalls that are watched by human operators, as well as activity logs make it difficult to send a great amount (MB, GB sized) data out.

### The Fake Al-Qaeda Pornography Arrest

![Maktab al-Dawa](http://www.biyokulule.com/sawiro/sawirada_waaweyn/Internet%20Jihad5_1.jpg)

Compare this to the arrest made on an **Al-Qaeda operative** who was found to have terrorist training manuals within  pornography videos in 2012. His method however, was less than sophisticated. It involved stuffing the TrueCrypt chunk within a forbidden zone of the video file, and altered the filesize of the video. Still, it is marvelous how much data was able to be fit within the file; megabyte/gigabyte.

The operative was "suspected" (right) to identify with a named "Maktab al-Dawa", a moniker not DIRECTLY affiliated with Al-Qaeda but used to spread political messages. Still, people argue that the "unnamed" assailant was just a lone-wolf sympathizer. This isn't as likely as a more likely fraudulent political scarecrow would likely do something like this; legally backed by the local authorities and covered by a false identity of witness protection agency. I MEAN COME ON, the messages claimed that "Terrorists are tired of law enforcement agencies being better than them, terrorists are tired of being poor, and that they may give up on plotting altogether thanks to successful force protection operations." LOL why would anyone put this on a porn video and shove it in his underwear? Who the hell is he even delivering it to?

Really, It's likely the guy was probably another *fed* trying to speak for the jihad jingoists and disappeared as soon as the arrest was made. I'm skeptical to this as well as many of the terrorist incidents in recent American media. Seems that there's a federal agency for talented actors, and you should really question the legitimacy of some of the more famous criminals and victims on the news. Most of these stories could be fabricated.

### So why should I use it?

Who knows?

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

**Recovered File in first frame after YouTube's H.264, no error code used.**
![](http://s24.postimg.org/vl6fn042d/final_result.jpg)

Filesize: 16KB, JPEG, BER = 000.29% Lena.jpg

**Recovered File in first frame hidden with error code r(50, 1) after YouTube**
![](http://s16.postimg.org/6myej8xjp/lena_image.jpg)

Filesize 16KB, JPEG, BER = 0% Lena.jpg

**Using error correcting codes makes embedding lossless so the file will be the exact same before and after. Luckily, the same parameters will work for all file/video combinations and there is no data dependency. PLEASE DO NOT EXPERIMENT WITH PARAMETERS UNLESS YOU KNOW WHAT YOU ARE DOING, YOU CAN CREATE DANGEROUS FILES. If you "brick" your computer running a corrupted program, I hate to say but I am not responsible.**

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
