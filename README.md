# ImagingCounting_RFID_InverseSolutions

The hardware system consists of several commercial off-the-shelf (COTS) radio-frequency identification (RFID) readers, and a multitude of passive RFID tags. The transmission between the reader and tags are described as follows:
1. The reader actively transmits modulated continuous-wave signals through passive antennas in the frequency range of 902-928 MHz. Baseband signal is modulated onto the carrier wave according to the Generation-2 (Gen2) Electronic Product Code (EPC) protocol. 50 channels, each with 500 kHz bandwidth, are successively used for transmission, where frequency hopping spread spectrum (FHSS) is adopted and the channel sequence is pre-determined by the random number generator in the reader.
2. The passive tags do not contain any battery, as they are able to harvest energy from the incoming reader's transmitted signal to power up the tag chip. The tag chip enables the tag to generate its own baseband signal (including tag ID and other command-specific information compliant with the Gen2 EPC protocol), and modulate its baseband signal onto the incoming reader's transmitted signal by changing its antenna's reflection coefficient. This process is called backscattering.
3. The reader receives the tag backscattering signal and decodes tag ID, as well as other RF information from the tags. Most current RFID readers support received signal strength indicator (RSSI) and phase information, which are the two most important RF parameters for this project.

This project uses RSSI and phase information available from the Impinj Speedway R420 reader, as well as tag ID, to construct the complex-valued reader-tag two-way transmission gain associated with combination of reader antenna, tag, and carrier frequency. Then the reflectivity image within the capture volume of interest is reconstructed using different inverse solutions.


This work has been published and please cite:
Guoyi Xu, Pragya Sharma, Xiaonan Hui, Edwin C. Kan, “Indoor object sensing using radio-frequency identification with inverse solutions”, 2020 IEEE Sensors Conference, Roterdam, Netherlands, Oct. 25-28, 2020.
