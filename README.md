GnuRadio-OCPC
=============

Oscillator Carrier phase compensation using Kalman tracking for  Gnu Radio  http://gnuradio.org

Contributors :  Venkhat.V (venkhat.ch28@gmail.com) , Kumar Karunakaran
In practice, Oscillator won’t produce carrier with accurate frequency fc, but within a small interval around fc. Hence, there is always a small frequency mismatch between carrier frequency from receiver’s local oscillator and transmitter carrier frequency. This causes the received symbol to get shifted by some phase (not magnitude). This phase have to be compensated before decoding.