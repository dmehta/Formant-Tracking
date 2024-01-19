karma.zip

Daryush D. Mehta, PhD
Center for Laryngeal Surgery and Voice Rehabilitation
Massachusetts General Hospital
Boston, MA
http://web.mit.edu/dmehta/www/

November 15, 2013

The purpose of this code is to estimate frequency and bandwidth values for formant and antiformant trajectories of a given speech waveform. A Kalman-based autoregressive moving average approach is employed [1-3].

The main routine is located in the 'run' folder: karma.m. Example usages are given in karma_demo.m with sample data located in the 'data' directory.

Please send questions or comments to Daryush Mehta (mehta.daryush@mgh.harvard.edu).

REFERENCES
----------

[1]   D. D. Mehta, D. Rudoy, and P. J. Wolfe, "Kalman-based autoregressive 
	moving average modeling and inference for formant and antiformant tracking," 
	The Journal of the Acoustical Society of America, vol. 132, no. 3, pp. 1732-1746, 2012.
[2]   D. Rudoy, D. N. Spendley, and P. J. Wolfe, "Conditionally linear
      Gaussian models for estimating vocal tract resonances," Proceedings of
      Interspeech, Antwerp, Belgium, 2007.
[3]   D. Rudoy, "Nonstationary time series modeling with application to speech
      signal processing," Doctor of Philosophy thesis, School of Engineering
      and Applied Sciences, Harvard University, Cambridge, MA, 2010.
      Chapter 3.

Version Info
------------

Version: 1.0, 04/10/2008 - Initial release
Version: 2.0, 07/31/2011 - Second release following manuscript submission
Version: 2.1, 11/15/2013 - Includes published manuscript reference (no code change from 2.0)

License
-------
The Formant Tracking software is ONLY available for academic purposes. It is not available for industrial or commercial applications of any kind without explicit arrangement with the author. The software must not be posted on any WWW or ftp sites or distributed in any other way without prior permission of the author. The author disclaims all warranties with regard to this software, including all implied warranties of merchantability and fitness. In no event shall the authors be liable for any special, indirect or consequential damages or any damages whatsoever resulting from loss of use, data or profits, whether in an action of contract, negligence or other tortious action, arising out of or in connection with the use or performance of this software. Permission to sell this software is not granted.

Installation
------------
Nothing fancy, just place in your favorite directory!