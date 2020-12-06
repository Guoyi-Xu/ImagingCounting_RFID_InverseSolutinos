package example;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.TimeoutException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.lang.Math;

import jxl.*;
import jxl.write.*;
import jxl.write.Number;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.jdom.JDOMException;

import org.llrp.ltk.exceptions.*;
import org.llrp.ltk.generated.enumerations.*;
import org.llrp.ltk.generated.interfaces.*;
import org.llrp.ltk.generated.messages.*;
import org.llrp.ltk.generated.parameters.*;
import org.llrp.ltk.types.*;
import org.llrp.ltk.net.LLRPConnection;
import org.llrp.ltk.net.LLRPConnectionAttemptFailedException;
import org.llrp.ltk.net.LLRPConnector;
import org.llrp.ltk.net.LLRPEndpoint;
import org.llrp.ltk.util.Util;

import org.llrp.ltk.generated.custom.enumerations.*;
import org.llrp.ltk.generated.custom.messages.*;
import org.llrp.ltk.generated.custom.parameters.*;

import com.jmatio.io.*;
import com.jmatio.common.*;
import com.jmatio.types.*;

public class MyDesignTestDp implements LLRPEndpoint {

	private LLRPConnection connection;

	private static Logger logger  = Logger.getLogger(MyDesignTestDp.class.getName());
  	private ROSpec rospec;
    private int MessageID = 15; // a random starting point
    private UnsignedInteger modelName;
    UnsignedShort maxPowerIndex;
    SignedShort maxPower;
    UnsignedShort maxSensIndex;
    SignedShort maxSens;
    UnsignedShort channelIndex;
    UnsignedShort hopTableID;
    static UnsignedShort antenna_sel;

    Double transmit_power;
    int receive_sensitivity;
    private static ArrayList<Double> transmitpower_table = new ArrayList<Double>();
    private static ArrayList<Integer> receivesensitivity_table = new ArrayList<Integer>();

    // create the parameters to be written into the excel
    private static ArrayList<Double> chindlist = new ArrayList<Double>();
    private static ArrayList<String> epclist = new ArrayList<String>();
    private static ArrayList<Double> antennalist = new ArrayList<Double>();
    private static ArrayList<Double> rssiimpinjlist = new ArrayList<Double>();
    private static ArrayList<Double> rssiimpinjlist_d = new ArrayList<Double>();
    private static ArrayList<Double> phasedeglist = new ArrayList<Double>();
    //private static ArrayList<Double> chindexlist = new ArrayList<Double>();
    private static ArrayList<Double> tagindexlist = new ArrayList<Double>();
    private static ArrayList<Double> msgfreqlist = new ArrayList<Double>();


    private static ArrayList<String> ftimelist = new ArrayList<String>();
    private static ArrayList<String> ltimelist = new ArrayList<String>();
    private static ArrayList<Double> ht = new ArrayList<Double>();
    private static ArrayList<Double> hti = new ArrayList<Double>();


    private static ArrayList<Double> nextlist = new ArrayList<Double>();

    private UnsignedInteger getUniqueMessageID() {
        return new UnsignedInteger(MessageID++);
    }

	public MyDesignTestDp() {
    }

	private void connect(String ip)
	{
		// create client-initiated LLRP connection
		connection = new LLRPConnector(this, ip);
		// connect to the reader
		// LLRPConnector.connect waits for successful
		// READER_EVENT_NOTIFICATION from reader
		try
		{
			logger.info("Initiate LLRP connection to reader");
			((LLRPConnector)connection).connect();
		}
		catch (LLRPConnectionAttemptFailedException e1)
		{
			e1.printStackTrace();
			System.exit(1);
		}
	}
	
	private void disconnect() {
		LLRPMessage response;
		CLOSE_CONNECTION close = new CLOSE_CONNECTION();
		close.setMessageID(getUniqueMessageID());
		try {
			// don't wait around too long for close
			response = connection.transact(close, 4000);
			
			// check whether the disconnection was successful
			StatusCode status = ((CLOSE_CONNECTION_RESPONSE)response).getLLRPStatus().getStatusCode();
			if (status.equals(new StatusCode("M_Success"))) {
				logger.info("CLOSE_CONNECTION was successful");
			}
			else {
				logger.info(response.toXMLString());
				logger.info("CLOSE_CONNECTION Failed ... continuing anyway");
			}
		} catch (InvalidLLRPMessageException ex) {
			logger.error("CLOSE_CONNECTION: Received invalid response message");
		} catch (TimeoutException ex) {
			logger.info("CLOSE_CONNECTION Timeouts ... continuing anyway");
		}
	}

	private void enableImpinjExtensions()
	{
		LLRPMessage response;
		try {
			logger.info("IMPINJ_ENABLE_EXTENSIONS ...");
			IMPINJ_ENABLE_EXTENSIONS ena = new IMPINJ_ENABLE_EXTENSIONS();
			ena.setMessageID(getUniqueMessageID());
			response = connection.transact(ena, 10000);
			StatusCode status = ((IMPINJ_ENABLE_EXTENSIONS_RESPONSE)response).getLLRPStatus().getStatusCode();
			if (status.equals(new StatusCode("M_Success"))) {
				logger.info("IMPINJ_ENABLE_EXTENSIONS was successful");
			}
			else {
				logger.info(response.toXMLString());
				logger.info("IMPINJ_ENABLE_EXTENSIONS Failure");
				System.exit(1);
			}
		}
		catch (InvalidLLRPMessageException ex) {
			logger.error("Could not process IMPINJ_ENABLE_EXTENSIONS response");
			System.exit(1);
		}
		catch (TimeoutException ex) {
			logger.error("Timeout Waiting for IMPINJ_ENABLE_EXTENSIONS response");
			System.exit(1);
		}
	}

	private void factoryDefault() {
		LLRPMessage response;
		try {
			// factory default the reader
			logger.info("SET_READER_CONFIG with factory default ...");
			SET_READER_CONFIG set = new SET_READER_CONFIG();
			set.setMessageID(getUniqueMessageID());
			set.setResetToFactoryDefault(new Bit(true));
			response = connection.transact(set, 10000);

			// check whether ROSpec addition was successful
			StatusCode status = ((SET_READER_CONFIG_RESPONSE)response).getLLRPStatus().getStatusCode();
			if (status.equals(new StatusCode("M_Success"))) {
				logger.info("SET_READER_CONFIG Factory Default was successful");
			}
			else {
				logger.info(response.toXMLString());
				logger.info("SET_READER_CONFIG Factory Default Failure");
				System.exit(1);
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void getReaderCapabilities() {
		LLRPMessage response;
		GET_READER_CAPABILITIES_RESPONSE gresp;
		GET_READER_CAPABILITIES get = new GET_READER_CAPABILITIES();
		GetReaderCapabilitiesRequestedData data = new GetReaderCapabilitiesRequestedData(GetReaderCapabilitiesRequestedData.All);
		get.setRequestedData(data);
		get.setMessageID(getUniqueMessageID());
		logger.info("Sending GET_READER_CAPABILITIES message ...");
		try {
			response = connection.transact(get, 10000);

			// check whether GET_CAPABILITIES addition was successful
			gresp = (GET_READER_CAPABILITIES_RESPONSE)response;
			StatusCode status = gresp.getLLRPStatus().getStatusCode();
			if (status.equals(new StatusCode("M_Success"))) {
				logger.info("GET_READER_CAPABILITIES was successful");

				// get the info we need
				GeneralDeviceCapabilities dev_cap = gresp.getGeneralDeviceCapabilities();
				if ((dev_cap == null) || (!dev_cap.getDeviceManufacturerName().equals(new UnsignedInteger(25882)))) {
					logger.error("MyDesignTestDp must use Impinj model Reader, not " + dev_cap.getDeviceManufacturerName().toString());
					System.exit(1);
				}

				modelName = dev_cap.getModelName();
				logger.info("Found Impinj reader model " + modelName.toString());

				// get the receive sensitivity table
				List<ReceiveSensitivityTableEntry> recv_list = dev_cap.getReceiveSensitivityTableEntryList();
				for (ReceiveSensitivityTableEntry recv_ent:recv_list) {
					String recv_ind = recv_ent.getIndex().toString();
					int recv_sens = recv_ent.getReceiveSensitivityValue().intValue();
					int recv_sens_real = recv_sens - 80;
					System.out.println(recv_ent.getIndex().intValue());
					receivesensitivity_table.add(recv_sens_real);
					System.out.println("Receive Sensitivity Index: " + recv_ind + ";");
					System.out.println("Receive Sensitivity Value: " + recv_sens_real + " dBm;");
				}

				// get the maximum receive sensitivity
				ReceiveSensitivityTableEntry recv_entry = recv_list.get(recv_list.size() - 1);
				maxSensIndex = recv_entry.getIndex();
				System.out.println("maxSensIndex is: " + maxSensIndex);
				maxSens = recv_entry.getReceiveSensitivityValue();
				System.out.println("max sens is: " + maxSens);
				System.out.println();

				// get the transmit power table
				if (gresp.getRegulatoryCapabilities() != null) {
					UHFBandCapabilities band_cap = gresp.getRegulatoryCapabilities().getUHFBandCapabilities();

					List<TransmitPowerLevelTableEntry> pwr_list = band_cap.getTransmitPowerLevelTableEntryList();
					for (TransmitPowerLevelTableEntry pwr:pwr_list) {
						String id = pwr.getIndex().toString();
						int p_v = pwr.getTransmitPowerValue().intValue();
						double p = p_v / 100.0;
						transmitpower_table.add(p);
						System.out.println("Power Index: " + id + "; Power Value: " + p + " dBm");
					}

					// Get the maximum transmit power (the last element in the transmit power level table).
					TransmitPowerLevelTableEntry entry = pwr_list.get(pwr_list.size() - 1);
					maxPowerIndex = entry.getIndex();
					maxPower = entry.getTransmitPowerValue();
					// LLRP sends power in dBm * 100;
					double d = ((double)maxPower.intValue())/100;
					logger.info("Max power " + d + " dBm at index " + maxPowerIndex.toString());

					// Print the frequency hopping table.
					FrequencyInformation fqinfo = band_cap.getFrequencyInformation();
					Bit hopping = fqinfo.getHopping();
					System.out.println("hopping or not: " + hopping);
					List<FrequencyHopTable> lhoptable = fqinfo.getFrequencyHopTableList();
					int sz_hoptable = lhoptable.size();
					System.out.println("Size of hoptable list: " + sz_hoptable);
					FrequencyHopTable hoptable = lhoptable.get(0);
					UnsignedIntegerArray fq_list = hoptable.getFrequency();
					int sz = fq_list.size();
					System.out.println("Size of frequency list: " + sz);
					int i = 0;
					while (i < fq_list.size()) {
						UnsignedInteger fq = fq_list.get(i);
						hti.add((i+1)/1.00);
						ht.add(fq.intValue()/1.00);
						System.out.println("frequency index: " + i + "; frequency value: " + fq);
						i += 1;
					}

				}

			}
			else {
				logger.info(response.toXMLString());
				logger.info("GET_READER_CAPABILITIES failures");
				System.exit(1);
			}
		} catch (InvalidLLRPMessageException ex) {
			logger.error("Could not display response string");
		} catch (TimeoutException ex) {
			logger.error("Timeout waiting for GET_READER_CAPABILITIES response");
			System.exit(1);
		}
	}

	private void getReaderConfiguration() {
		LLRPMessage response;
		GET_READER_CONFIG_RESPONSE gresp;
		GET_READER_CONFIG get = new GET_READER_CONFIG();
		GetReaderConfigRequestedData data = new GetReaderConfigRequestedData(GetReaderConfigRequestedData.All);
		get.setRequestedData(data);
		get.setMessageID(getUniqueMessageID());
		get.setAntennaID(new UnsignedShort(0));
		get.setGPIPortNum(new UnsignedShort(0));
		get.setGPOPortNum(new UnsignedShort(0));

		logger.info("Sending GET_READER_CONFIG message ...");
		try {
			response = connection.transact(get, 10000);

			// check whether GET_CONFIG addition was successful
			gresp = (GET_READER_CONFIG_RESPONSE)response;
			StatusCode status = gresp.getLLRPStatus().getStatusCode();
			if (status.equals(new StatusCode("M_Success"))) {
				logger.info("GET_READER_CONFIG was successful");

				// print the antenna configurations
				List<AntennaConfiguration> alist = gresp.getAntennaConfigurationList();
				System.out.println(alist.size());
				if (!alist.isEmpty()) {
					for (int j = 0; j < alist.size(); j ++) {
						AntennaConfiguration a_cfg = alist.get(j);
						UnsignedShort aid = a_cfg.getAntennaID();
						channelIndex = a_cfg.getRFTransmitter().getChannelIndex();
						hopTableID = a_cfg.getRFTransmitter().getHopTableID();
						System.out.println("Antenna ID is: " + aid);
						logger.info("ChannelIndex " + channelIndex.toString() + "\nhopTableID: " + hopTableID.toString());

						// Print the transmit power and receive sensitivity configuration
						UnsignedShort transmit_power_index = a_cfg.getRFTransmitter().getTransmitPower();
						transmit_power = transmitpower_table.get(transmit_power_index.intValue() - 1);
						UnsignedShort receive_sensitivity_index = a_cfg.getRFReceiver().getReceiverSensitivity();
						receive_sensitivity = receivesensitivity_table.get(receive_sensitivity_index.intValue() - 1);
						System.out.println("Transmit Power Table index is: " + transmit_power_index + "\nTransmit Power is: " + transmit_power + " dBm;");
						System.out.println("Receive Sensitivity Table index is: " + receive_sensitivity_index + "\nReceive Sensitivity is: " + receive_sensitivity + " dBm;");

					}				
				} else {
					logger.error("Could not find antenna configuration.");
					System.exit(1);
				}

				// print the antenna properties
				List<AntennaProperties> plist = gresp.getAntennaPropertiesList();
				System.out.println(plist.size());
				if (!plist.isEmpty()) {
					for (int i = 0; i < plist.size(); i ++) {
						AntennaProperties a_pro = plist.get(i);
						UnsignedShort a_id = a_pro.getAntennaID();
						SignedShort a_gain = a_pro.getAntennaGain();
						Bit check = a_pro.getAntennaConnected();
						System.out.println("Antenna ID is: " + a_id.toString());
						System.out.println("Antenna gain is: " + a_gain.toString());
						System.out.println("The antenna is connected? " + check);
					}
				} else {
					logger.error("Could not find antenna properties.");
					System.exit(1);
				}


			}
			else {
				logger.info(response.toXMLString());
				logger.info("GET_READER_CONFIG failures");
				System.exit(1);
			}
		} catch (InvalidLLRPMessageException ex) {
			logger.error("Could not display response string");
		} catch (TimeoutException ex) {
			logger.error("Timeout waiting for GET_READER_CONFIG response");
			System.exit(1);
		}
	}

	private void setReaderConfiguration() {
		LLRPMessage response;

		logger.info("Loading SET_READER_CONFIG message from file SET_READER_CONFIG.xml ...");
		try {
			LLRPMessage setConfigMsg = Util.loadXMLLLRPMessage(new File("./source/example/SET_READER_CONFIG.xml"));
			// TODO make sure this is a SET_READER_CONFIG message
			
			// touch up the transmit power for max
			SET_READER_CONFIG setConfig = (SET_READER_CONFIG)setConfigMsg;
			AntennaConfiguration a_cfg = setConfig.getAntennaConfigurationList().get(0);
			RFTransmitter rftx = a_cfg.getRFTransmitter();
			rftx.setChannelIndex(channelIndex);
			rftx.setHopTableID(hopTableID);
			rftx.setTransmitPower(maxPowerIndex);

			// set the receive sensitivity to max
			RFReceiver rfrx = a_cfg.getRFReceiver();
			//rfrx.setReceiverSensitivity(maxSensIndex);

			// set the reader event notification spec
			ReaderEventNotificationSpec rens = new ReaderEventNotificationSpec();
			List<EventNotificationState> enslist = new ArrayList<EventNotificationState>();
			EventNotificationState ens = new EventNotificationState();
			ens.setNotificationState(new Bit(1));
			ens.setEventType(new NotificationEventType(NotificationEventType.Upon_Hopping_To_Next_Channel));
			enslist.add(ens);
			rens.setEventNotificationStateList(enslist);
			setConfig.setReaderEventNotificationSpec(rens);
			
			response = connection.transact(setConfig, 10000);
			
			// check whether SET_READER_CONFIG addition was successful
			StatusCode status = ((SET_READER_CONFIG_RESPONSE)response).getLLRPStatus().getStatusCode();
			if (status.equals(new StatusCode("M_Success"))) {
				logger.info("SET_READER_CONFIG was successful");
			}
			else {
				logger.info(response.toXMLString());
				logger.info("SET_READER_CONFIG failures");
				System.exit(1);
			}
		} catch (TimeoutException ex) {
			logger.error("Timeout waiting for SET_READER_CONFIG response");
			System.exit(1);
		} catch (FileNotFoundException ex) {
			logger.error("Could not find file");
			System.exit(1);
		} catch (IOException ex) {
			logger.error("IO Exception on file");
			System.exit(1);
		} catch (JDOMException ex) {
			logger.error("Unable to convert LTK-XML to DOM");
			System.exit(1);
		} catch (InvalidLLRPMessageException ex) {
			logger.error("Unable to convert LTK-XML to Internal Object");
			System.exit(1);
		}
	}

	private ADD_ROSPEC buildROSpecFromObjects(UnsignedShort ant_sel) {
		logger.info("Building ADD_ROSPEC message from scratch ...");
		ADD_ROSPEC addRoSpec = new ADD_ROSPEC();
		addRoSpec.setMessageID(getUniqueMessageID());
		rospec = new ROSpec();

		// set up the basic info for the RO Spec.
		rospec.setCurrentState(new ROSpecState(ROSpecState.Disabled));
		rospec.setPriority(new UnsignedByte(0));
		rospec.setROSpecID(new UnsignedInteger(1));

		// set the start and stop conditions for the ROSpec.
		// For now, we will start and stop manually.
		ROBoundarySpec boundary = new ROBoundarySpec();
		ROSpecStartTrigger start = new ROSpecStartTrigger();
		ROSpecStopTrigger stop = new ROSpecStopTrigger();
		start.setROSpecStartTriggerType(new ROSpecStartTriggerType(ROSpecStartTriggerType.Null));
		stop.setROSpecStopTriggerType(new ROSpecStopTriggerType(ROSpecStopTriggerType.Null));
		stop.setDurationTriggerValue(new UnsignedInteger(0));
		boundary.setROSpecStartTrigger(start);
		boundary.setROSpecStopTrigger(stop);
		rospec.setROBoundarySpec(boundary);

		// set up what we want to do in the ROSpec. In this case
		// build the simples inventory on all channels using defaults
		AISpec aispec = new AISpec();
		// what antennas to use.
		UnsignedShortArray ants = new UnsignedShortArray();
		ants.add(ant_sel);	// 0 means all antennas
		aispec.setAntennaIDs(ants);
		// set up the AISpec stop condition and options for inventory
		AISpecStopTrigger aistop = new AISpecStopTrigger();
		aistop.setAISpecStopTriggerType(new AISpecStopTriggerType(AISpecStopTriggerType.Null));
		aistop.setDurationTrigger(new UnsignedInteger(0));
		aispec.setAISpecStopTrigger(aistop);
		// set up any override configuration. none in this case
		InventoryParameterSpec ispec = new InventoryParameterSpec();
		ispec.setAntennaConfigurationList(null);
		ispec.setInventoryParameterSpecID(new UnsignedShort(1));
		ispec.setProtocolID(new AirProtocols(AirProtocols.EPCGlobalClass1Gen2));
		List<InventoryParameterSpec> ilist = new ArrayList<InventoryParameterSpec>();
		ilist.add(ispec);
		aispec.setInventoryParameterSpecList(ilist);
		List<SpecParameter> slist = new ArrayList<SpecParameter>();
		slist.add(aispec);

		rospec.setSpecParameterList(slist);

		// Set up when we want to see the report.
		ROReportSpec report = new ROReportSpec();
		report.setROReportTrigger(new ROReportTriggerType(ROReportTriggerType.Upon_N_Tags_Or_End_Of_ROSpec));
		report.setN(new UnsignedShort(1));
		// Set up the LLRP tag report configurations.
		TagReportContentSelector tagrcs = new TagReportContentSelector();
		tagrcs.setEnableROSpecID(new Bit(false));
		tagrcs.setEnableSpecIndex(new Bit(false));
		tagrcs.setEnableInventoryParameterSpecID(new Bit(false));
		tagrcs.setEnableAntennaID(new Bit(true));
		tagrcs.setEnableChannelIndex(new Bit(true));
		tagrcs.setEnablePeakRSSI(new Bit(false));
		tagrcs.setEnableFirstSeenTimestamp(new Bit(true));
		tagrcs.setEnableLastSeenTimestamp(new Bit(true));
		tagrcs.setEnableTagSeenCount(new Bit(false));
		tagrcs.setEnableAccessSpecID(new Bit(false));
		report.setTagReportContentSelector(tagrcs);

		// Set up the Impinj tag report configurations.
		List<Custom> itagrlist = new ArrayList<Custom>();
		ImpinjTagReportContentSelector itagrcs = new ImpinjTagReportContentSelector();
		// Set the Peak RSSI Enabled.
		ImpinjEnablePeakRSSI rssi = new ImpinjEnablePeakRSSI();
		rssi.setPeakRSSIMode(new ImpinjPeakRSSIMode(ImpinjPeakRSSIMode.Enabled));
		itagrcs.setImpinjEnablePeakRSSI(rssi);
		// Set the RF Phase Angle Enabled.
		ImpinjEnableRFPhaseAngle phase = new ImpinjEnableRFPhaseAngle();
		phase.setRFPhaseAngleMode(new ImpinjRFPhaseAngleMode(ImpinjRFPhaseAngleMode.Enabled));
		itagrcs.setImpinjEnableRFPhaseAngle(phase);
		// Set the TID Enabled.
		ImpinjEnableSerializedTID tid = new ImpinjEnableSerializedTID();
		tid.setSerializedTIDMode(new ImpinjSerializedTIDMode(ImpinjSerializedTIDMode.Disabled));
		itagrcs.setImpinjEnableSerializedTID(tid);
		itagrlist.add((Custom)itagrcs);
		report.setCustomList(itagrlist);

		rospec.setROReportSpec(report);
		addRoSpec.setROSpec(rospec);
		return addRoSpec;
	}

	private ADD_ROSPEC buildROSpecFromFile() {
		logger.info("Loading ADD_ROSPEC message from file ADD_ROSPEC.xml ...");
		try {
			LLRPMessage addRospec = Util.loadXMLLLRPMessage(new File("./source/example/ADD_ROSPEC.xml"));
			return (ADD_ROSPEC) addRospec;
		} catch (FileNotFoundException ex) {
			logger.error("Could not find file");
			System.exit(1);
		} catch (IOException ex) {
			logger.error("IO Exception on file");
			System.exit(1);
		} catch (JDOMException ex) {
			logger.error("Unable to convert LTK-XML to DOM");
			System.exit(1);
		} catch (InvalidLLRPMessageException ex) {
			logger.error("Unable to convert LTK-XML to Internal Object");
			System.exit(1);
		}
		return null;
	}

	private void addRoSpec(boolean xml, UnsignedShort ant_sel) {
        LLRPMessage response;

        ADD_ROSPEC addRospec = null;

        if(xml) {
            addRospec = buildROSpecFromFile();
        } else {
            addRospec = buildROSpecFromObjects(ant_sel);
        }

        addRospec.setMessageID(getUniqueMessageID());
        rospec = addRospec.getROSpec();
        
        logger.info("Sending ADD_ROSPEC message  ...");
        try {
            response =  connection.transact(addRospec, 10000);

            // check whether ROSpec addition was successful
            StatusCode status = ((ADD_ROSPEC_RESPONSE)response).getLLRPStatus().getStatusCode();
            if (status.equals(new StatusCode("M_Success"))) {
                    logger.info("ADD_ROSPEC was successful");
            }
            else {
                    logger.info(response.toXMLString());
                    logger.info("ADD_ROSPEC failures");
                    System.exit(1);
            }
        } catch (InvalidLLRPMessageException ex) {
            logger.error("Could not display response string");
        } catch (TimeoutException ex) {
            logger.error("Timeout waiting for ADD_ROSPEC response");
            System.exit(1);
        }
    }

	private void enable() {
		LLRPMessage response;
		try {
			// factory default the reader
			logger.info("ENABLE_ROSPEC ...");
			ENABLE_ROSPEC ena = new ENABLE_ROSPEC();
			ena.setMessageID(getUniqueMessageID());
			ena.setROSpecID(rospec.getROSpecID());
			response = connection.transact(ena, 10000);
			// check whether ROSpec addition was successful
			StatusCode status = ((ENABLE_ROSPEC_RESPONSE) response).getLLRPStatus().getStatusCode();
			if (status.equals(new StatusCode("M_Success"))) {
				logger.info("ENABLE_ROSPEC was successful");
			}
			else {
				logger.error(response.toXMLString());
				logger.info("ENABLE_ROSPEC_RESPONSE failed");
				System.exit(1);
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void start() {
		LLRPMessage response;
		try {
			logger.info("START_ROSPEC ...");
			START_ROSPEC start = new START_ROSPEC();
			start.setMessageID(getUniqueMessageID());
			start.setROSpecID(rospec.getROSpecID());
			response = connection.transact(start, 10000);
			// check whether ROSpec addition was successful
			StatusCode status = ((START_ROSPEC_RESPONSE) response).getLLRPStatus().getStatusCode();
			if (status.equals(new StatusCode("M_Success"))) {
				logger.info("START_ROSPEC was successful");
			}
			else {
				logger.error(response.toXMLString());
				logger.info("START_ROSPEC_RESPONSE failed");
				System.exit(1);
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void stop() {
		LLRPMessage response;
		try {
			logger.info("STOP_ROSPEC ...");
			STOP_ROSPEC stop = new STOP_ROSPEC();
			stop.setMessageID(getUniqueMessageID());
			stop.setROSpecID(rospec.getROSpecID());

			response = connection.transact(stop, 10000);

			// check whether ROSpec addition was successful
			StatusCode status = ((STOP_ROSPEC_RESPONSE)response).getLLRPStatus().getStatusCode();
			if (status.equals(new StatusCode("M_Success"))) {
				logger.info("STOP_ROSPEC was successful");
			} else {
				logger.error(response.toXMLString());
				logger.info("STOP_ROSPEC_RESPONSE failed");
				System.exit(1);
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void delete() {
		LLRPMessage response;
		try {
			logger.info("DELETE_ROSPEC ...");
			DELETE_ROSPEC delete = new DELETE_ROSPEC();
			delete.setMessageID(getUniqueMessageID());
			delete.setROSpecID(rospec.getROSpecID());

			response = connection.transact(delete, 10000);

			// check whether ROSpec addition was successful
			StatusCode status = ((DELETE_ROSPEC_RESPONSE)response).getLLRPStatus().getStatusCode();
			if (status.equals(new StatusCode("M_Success"))) {
				logger.info("DELETE_ROSPEC was successful");
			} else {
				logger.error(response.toXMLString());
				logger.info("DELETE_ROSPEC_RESPONSE failed");
				System.exit(1);
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void logOneCustom(Custom cust) {

        String output = "";
        if(!cust.getVendorIdentifier().equals(25882)) {
            logger.error("Non Impinj Extension Found in message");
            return;
        }
        LLRPParameter report_cust = (LLRPParameter)cust;
    }


	private void logOneTagReport(TagReportData tr) {
		// As an example here, we will just get the stuff out of here and
		// for a super long string
		LLRPParameter epcp = (LLRPParameter) tr.getEPCParameter();
		// epc is not optional, so we should fail if we can't find it
		String epcString = "EPC: ";
		if (epcp != null) {
			if (epcp.getName().equals("EPC_96")) {
				EPC_96 epc96 = (EPC_96) epcp;
				epcString += epc96.getEPC().toString() + "\n";
			} else if (epcp.getName().equals("EPCData")) {
				EPCData epcData = (EPCData) epcp;
				epcString += epcData.getEPC().toString() + "\n";
			}
		} else {
			logger.error("Could not find EPC in Tag Report");
			System.exit(1);
		}
		// all of these values are optional, so check their non-nullness first
		String antenna = "unknownAntenna";
		if (tr.getAntennaID() != null) {
			antenna = tr.getAntennaID().getAntennaID().toString();
			epcString += "Antenna: " + antenna + "\n";
		}
		String index = "unknownChannelIndex";
		if (tr.getChannelIndex() != null) {
			index = tr.getChannelIndex().getChannelIndex().toString();
			epcString += "ChanIndex: " + index + "\n";
		}
		String ftime = "unknownFirstTimestamp";
		if (tr.getFirstSeenTimestampUTC() != null) {
			ftime = tr.getFirstSeenTimestampUTC().getMicroseconds().toString();
			epcString += "FirstSeen: " + ftime + "\n";
		}
		String paramId = "unknownParamSpecID";
		if (tr.getInventoryParameterSpecID() != null) {
			paramId = tr.getInventoryParameterSpecID().getInventoryParameterSpecID().toString();
			epcString += "ParamSpecID: " + paramId + "\n";
		}
		String ltime = "unknownLastTimestamp";
		if (tr.getLastSeenTimestampUTC() != null) {
			ltime = tr.getLastSeenTimestampUTC().getMicroseconds().toString();
			epcString += "LastSeen: " + ltime + "\n";
		}
		String rssi = "unknownRssi";
		if (tr.getPeakRSSI() != null) {
			rssi = tr.getPeakRSSI().getPeakRSSI().toString();
			epcString += "RSSI: " + rssi + "\n";
		}
		String roid = "unknownRospecID";
		if (tr.getROSpecID() != null) {
			roid = tr.getROSpecID().getROSpecID().toString();
			epcString += "ROSpecID: " + roid + "\n";
		}
		String tagseen = "unknownTagSeenCount";
		if (tr.getTagSeenCount() != null) {
			tagseen = tr.getTagSeenCount().getTagCount().toString();
			epcString += "SeenCount: " + tagseen + "\n";
		}
		String ftimeup = "unknownFirstTimeUptime";
		if (tr.getTagSeenCount() != null) {
			tagseen = tr.getFirstSeenTimestampUptime().getMicroseconds().toString();
			epcString += "SeenCount: " + tagseen + "\n";
		}
		String ltimeup = "unknownLastTimeUptime";
		if (tr.getTagSeenCount() != null) {
			tagseen = tr.getLastSeenTimestampUptime().getMicroseconds().toString();
			epcString += "SeenCount: " + tagseen + "\n";
		}
		
		System.out.print(epcString);

		List<Custom> clist = tr.getCustomList();
		Double rssiimpinj;
		Double rssiimpinj_d;
		int phase;
		Double phase_deg;

		Double msg_freq = ht.get(tr.getChannelIndex().getChannelIndex().intValue() - 1);
		msgfreqlist.add(msg_freq);

		chindlist.add(tr.getChannelIndex().getChannelIndex().toInteger()/1.00);
		if (epcp.getName().equals("EPC_96")) {
			EPC_96 epc96 = (EPC_96) epcp;
			epclist.add(epc96.getEPC().toString());
		} else if (epcp.getName().equals("EPCData")) {
			EPCData epcData = (EPCData) epcp;
			epclist.add(epcData.getEPC().toString());
		}
		String epc_str = epclist.get(epclist.size() - 1);
		if (epc_str.contains("300833B2DDD9014000004001")) {
			tagindexlist.add(new Double(1));
		} else if (epc_str.contains("300833B2DDD9014000004002")) {
			tagindexlist.add(new Double(2));
		} else if (epc_str.contains("300833B2DDD9014000004003")) {
			tagindexlist.add(new Double(3));
		} else if (epc_str.contains("300833B2DDD9014000004004")) {
			tagindexlist.add(new Double(4));
		} else if (epc_str.contains("300833B2DDD9014000004005")) {
			tagindexlist.add(new Double(5));
		} else if (epc_str.contains("300833B2DDD9014000004006")) {
			tagindexlist.add(new Double(6));
		} else if (epc_str.contains("300833B2DDD9014000004007")) {
			tagindexlist.add(new Double(7));
		} else if (epc_str.contains("300833B2DDD9014000004008")) {
			tagindexlist.add(new Double(8));
		} else if (epc_str.contains("300833B2DDD9014000004009")) {
			tagindexlist.add(new Double(9));
		} else if (epc_str.contains("300833B2DDD9014000004010")) {
			tagindexlist.add(new Double(10));
		} else if (epc_str.contains("300833B2DDD9014000004011")) {
			tagindexlist.add(new Double(11));
		} else if (epc_str.contains("300833B2DDD9014000004012")) {
			tagindexlist.add(new Double(12));
		} else if (epc_str.contains("300833B2DDD9014000004013")) {
			tagindexlist.add(new Double(13));
		} else if (epc_str.contains("300833B2DDD9014000004014")) {
			tagindexlist.add(new Double(14));
		} else if (epc_str.contains("300833B2DDD9014000004015")) {
			tagindexlist.add(new Double(15));
		} else if (epc_str.contains("300833B2DDD9014000004016")) {
			tagindexlist.add(new Double(16));
		} else if (epc_str.contains("300833B2DDD9014000004017")) {
			tagindexlist.add(new Double(17));
		} else if (epc_str.contains("300833B2DDD9014000004018")) {
			tagindexlist.add(new Double(18));
		} else if (epc_str.contains("300833B2DDD9014000004019")) {
			tagindexlist.add(new Double(19));
		} else if (epc_str.contains("300833B2DDD9014000004020")) {
			tagindexlist.add(new Double(20));
		} else if (epc_str.contains("300833B2DDD9014000004021")) {
			tagindexlist.add(new Double(21));
		} else if (epc_str.contains("300833B2DDD9014000004022")) {
			tagindexlist.add(new Double(22));
		} else if (epc_str.contains("300833B2DDD9014000004023")) {
			tagindexlist.add(new Double(23));
		} else if (epc_str.contains("300833B2DDD9014000004024")) {
			tagindexlist.add(new Double(24));
		} else if (epc_str.contains("300833B2DDD9014000004025")) {
			tagindexlist.add(new Double(25));
		} else if (epc_str.contains("300833B2DDD9014000004026")) {
			tagindexlist.add(new Double(26));
		} else if (epc_str.contains("300833B2DDD9014000004027")) {
			tagindexlist.add(new Double(27));
		} else if (epc_str.contains("300833B2DDD9014000004028")) {
			tagindexlist.add(new Double(28));
		} else if (epc_str.contains("300833B2DDD9014000004029")) {
			tagindexlist.add(new Double(29));
		} else if (epc_str.contains("300833B2DDD9014000004030")) {
			tagindexlist.add(new Double(30));
		} else if (epc_str.contains("300833B2DDD9014000004031")) {
			tagindexlist.add(new Double(31));
		} else if (epc_str.contains("300833B2DDD9014000004032")) {
			tagindexlist.add(new Double(32));
		} else if (epc_str.contains("300833B2DDD9014000004033")) {
			tagindexlist.add(new Double(33));
		} else if (epc_str.contains("300833B2DDD9014000004034")) {
			tagindexlist.add(new Double(34));
		} else if (epc_str.contains("300833B2DDD9014000004035")) {
			tagindexlist.add(new Double(35));
		} else if (epc_str.contains("300833B2DDD9014000004036")) {
			tagindexlist.add(new Double(36));
		} else if (epc_str.contains("300833B2DDD9014000004037")) {
			tagindexlist.add(new Double(37));
		} else if (epc_str.contains("300833B2DDD9014000004038")) {
			tagindexlist.add(new Double(38));
		} else if (epc_str.contains("300833B2DDD9014000004039")) {
			tagindexlist.add(new Double(39));
		} else if (epc_str.contains("300833B2DDD9014000004040")) {
			tagindexlist.add(new Double(40));
		} else if (epc_str.contains("300833B2DDD9014000004041")) {
			tagindexlist.add(new Double(41));
		} else if (epc_str.contains("300833B2DDD9014000004042")) {
			tagindexlist.add(new Double(42));
		} else if (epc_str.contains("300833B2DDD9014000004043")) {
			tagindexlist.add(new Double(43));
		} else if (epc_str.contains("300833B2DDD9014000004044")) {
			tagindexlist.add(new Double(44));
		} else if (epc_str.contains("300833B2DDD9014000004045")) {
			tagindexlist.add(new Double(45));
		} else if (epc_str.contains("300833B2DDD9014000004046")) {
			tagindexlist.add(new Double(46));
		} else if (epc_str.contains("300833B2DDD9014000004047")) {
			tagindexlist.add(new Double(47));
		} else if (epc_str.contains("300833B2DDD9014000004048")) {
			tagindexlist.add(new Double(48));
		} else if (epc_str.contains("300833B2DDD9014000004049")) {
			tagindexlist.add(new Double(49));
		} else if (epc_str.contains("300833B2DDD9014000004050")) {
			tagindexlist.add(new Double(50));
		} else if (epc_str.contains("300833B2DDD9014000004051")) {
			tagindexlist.add(new Double(51));
		} else if (epc_str.contains("300833B2DDD9014000004052")) {
			tagindexlist.add(new Double(52));
		} else if (epc_str.contains("300833B2DDD9014000004053")) {
			tagindexlist.add(new Double(53));
		} else if (epc_str.contains("300833B2DDD9014000004054")) {
			tagindexlist.add(new Double(54));
		} else if (epc_str.contains("300833B2DDD9014000004055")) {
			tagindexlist.add(new Double(55));
		} else if (epc_str.contains("300833B2DDD9014000004056")) {
			tagindexlist.add(new Double(56));
		} else if (epc_str.contains("300833B2DDD9014000004057")) {
			tagindexlist.add(new Double(57));
		} else if (epc_str.contains("300833B2DDD9014000004058")) {
			tagindexlist.add(new Double(58));
		} else if (epc_str.contains("300833B2DDD9014000004059")) {
			tagindexlist.add(new Double(59));
		} else if (epc_str.contains("300833B2DDD9014000004060")) {
			tagindexlist.add(new Double(60));
		} else if (epc_str.contains("300833B2DDD9014000004061")) {
			tagindexlist.add(new Double(61));
		} else if (epc_str.contains("300833B2DDD9014000004062")) {
			tagindexlist.add(new Double(62));
		} else if (epc_str.contains("300833B2DDD9014000004063")) {
			tagindexlist.add(new Double(63));
		} else if (epc_str.contains("300833B2DDD9014000004064")) {
			tagindexlist.add(new Double(64));
		} else if (epc_str.contains("300833B2DDD9014000004065")) {
			tagindexlist.add(new Double(65));
		} else if (epc_str.contains("300833B2DDD9014000004066")) {
			tagindexlist.add(new Double(66));
		} else if (epc_str.contains("300833B2DDD9014000004067")) {
			tagindexlist.add(new Double(67));
		} else if (epc_str.contains("300833B2DDD9014000004068")) {
			tagindexlist.add(new Double(68));
		} else if (epc_str.contains("300833B2DDD9014000004069")) {
			tagindexlist.add(new Double(69));
		} else if (epc_str.contains("300833B2DDD9014000004070")) {
			tagindexlist.add(new Double(70));
		} else if (epc_str.contains("300833B2DDD9014000004071")) {
			tagindexlist.add(new Double(71));
		} else if (epc_str.contains("300833B2DDD9014000004072")) {
			tagindexlist.add(new Double(72));
		} else if (epc_str.contains("300833B2DDD9014000004073")) {
			tagindexlist.add(new Double(73));
		} else if (epc_str.contains("300833B2DDD9014000004074")) {
			tagindexlist.add(new Double(74));
		} else if (epc_str.contains("300833B2DDD9014000004075")) {
			tagindexlist.add(new Double(75));
		} else if (epc_str.contains("300833B2DDD9014000004076")) {
			tagindexlist.add(new Double(76));
		} else if (epc_str.contains("300833B2DDD9014000004077")) {
			tagindexlist.add(new Double(77));
		} else if (epc_str.contains("300833B2DDD9014000004078")) {
			tagindexlist.add(new Double(78));
		} else if (epc_str.contains("300833B2DDD9014000004079")) {
			tagindexlist.add(new Double(79));
		} else if (epc_str.contains("300833B2DDD9014000004080")) {
			tagindexlist.add(new Double(80));
		} else {
			tagindexlist.add(new Double(0));
		}
		antennalist.add(tr.getAntennaID().getAntennaID().toInteger()/1.00);
		for (Custom cd: clist) {
			if (cd.getClass() == ImpinjPeakRSSI.class) {
				rssiimpinj = ((ImpinjPeakRSSI)cd).getRSSI().toInteger()/100.000;
				rssiimpinjlist.add(rssiimpinj);
				rssiimpinj_d = Math.pow(10, rssiimpinj/10)*1000000;
				rssiimpinjlist_d.add(rssiimpinj_d);
				System.out.println("RSSIimpinj: " + rssiimpinj);
				System.out.println("RSSIimpinj: " + rssiimpinj_d);
			}
			if (cd.getClass() == ImpinjRFPhaseAngle.class) {
				phase = ((ImpinjRFPhaseAngle)cd).getPhaseAngle().toInteger();
				phase_deg = 360.0*phase/4096.0;
				phasedeglist.add(phase_deg);
				System.out.println("RFPhase: " + phase_deg);
			}
		}
		ftimelist.add(tr.getFirstSeenTimestampUTC().getMicroseconds().toString());
		ltimelist.add(tr.getLastSeenTimestampUTC().getMicroseconds().toString());

		System.out.println();
	}


	private void logReaderEvent(ReaderEventNotificationData rend) {
		if (rend.getHoppingEvent() != null) {
			HoppingEvent hp = rend.getHoppingEvent();
			UnsignedShort hoptableID = hp.getHopTableID();
			int next = hp.getNextChannelIndex().toInteger().intValue();
			nextlist.add((Double)(next/1.00));
			System.out.println("Hop Table ID is: " + hoptableID);
			System.out.println("Next Channel Index is: " + next);
			System.out.println();
		} else {
			logger.error("Could not find hopping event. ");
		}
	}

	// messageReceived method is called whenever a message is received
	// asynchronously on the LLRP connection.
	public void messageReceived(LLRPMessage message) {
		// convert all messages received to LTK-XML representation
		// and print them to the console
		logger.debug("Received" + message.getName() + " message asynchronously");
		System.out.println("Message received.");

		if (message.getTypeNum() == RO_ACCESS_REPORT.TYPENUM) {
			RO_ACCESS_REPORT report = (RO_ACCESS_REPORT) message;
			List<TagReportData> tdlist = report.getTagReportDataList();
			for (TagReportData tr : tdlist) {
				logOneTagReport(tr);
			}

			List<Custom> clist = report.getCustomList();
            for (Custom cust : clist) {
                logOneCustom(cust);
            }

		} else if (message.getTypeNum() == READER_EVENT_NOTIFICATION.TYPENUM) {
			READER_EVENT_NOTIFICATION report = (READER_EVENT_NOTIFICATION) message;
			ReaderEventNotificationData rend = report.getReaderEventNotificationData();
			logReaderEvent(rend);
		}
	}
	
	@Override
	public void errorOccured(String arg0) {
		// TODO Auto-generated method stub
		
	}	

    /**
     * @param args
     */
    public static void main(String[] args) throws IOException, WriteException
    {
    	BasicConfigurator.configure();

    	if (args.length < 1) {
    		System.out.print("Must pass reader hostname or IP as argument 1");
    		System.exit(-1);
    	}

    	Logger.getRootLogger().setLevel(Level.ERROR);
    	MyDesignTestDp masterpiece = new MyDesignTestDp();
    	logger.setLevel(Level.INFO);

    	String usecase = new String("Center");
    	String caseindex = new String("1");
    	// Part I: Calibration.
    	masterpiece.connect(args[0]);
    	masterpiece.enableImpinjExtensions();
    	masterpiece.factoryDefault();
        masterpiece.getReaderCapabilities();
        masterpiece.getReaderConfiguration();
        masterpiece.setReaderConfiguration();
        //masterpiece.getReaderConfiguration();

        // start of for loop
    	int cnt = 1;
    	while (cnt <= 1) {
    		// Record the collected data variable size before the next update period.
    		int s1 = chindlist.size();

    		// Switch among antennas.
    		for (int i = 3; i < 4; i ++) {
        		switch (i) {
        			case 1: antenna_sel = new UnsignedShort(1); break;
        			case 2: antenna_sel = new UnsignedShort(2); break;
        			case 3: antenna_sel = new UnsignedShort(3); break;
        			case 4: antenna_sel = new UnsignedShort(4); break;
        			default: break;
        		}
        		masterpiece.addRoSpec(false, antenna_sel);
    			masterpiece.enable();
    			masterpiece.start();

    			try {
    				Thread.sleep(60000);
    			} catch (InterruptedException ex) {
    				logger.error("Sleep Interruped");
    			}

    			masterpiece.stop();
    			masterpiece.delete();
    			
        	}

        	// Record the collected data variable size after running the current update period.
        	int s2 = chindlist.size();
    		int count = s2 - s1;
        	//int count = chindexlist.size();
    		double[] clist_d = new double[count];
    		double[] tlist_d = new double[count];
    		double[] alist_d = new double[count];
    		double[] rlist_d = new double[count];
    		double[] rdlist_d = new double[count];
    		double[] plist_d = new double[count];
    		double[] flist_d = new double[count];

    		System.out.println("s1 is: " + s1);
    		System.out.println("s2 is: " + s2);

    		for (int i = s1; i < s2; i ++) {
    			clist_d[i - s1] = chindlist.get(i);
    			tlist_d[i - s1] = tagindexlist.get(i);
    			alist_d[i - s1] = antennalist.get(i);
    			rlist_d[i - s1] = rssiimpinjlist.get(i);
    			rdlist_d[i - s1] = rssiimpinjlist_d.get(i);
    			plist_d[i - s1] = phasedeglist.get(i);
    			flist_d[i - s1] = msgfreqlist.get(i);
    		}
    	
    		MLDouble mlDouble1 = new MLDouble("chindlist", clist_d, count);
    		MLDouble mlDouble2 = new MLDouble("tagindexlist", tlist_d, count);
    		MLDouble mlDouble3 = new MLDouble("antennalist", alist_d, count);
    		MLDouble mlDouble4 = new MLDouble("rssiimpinjlist", rlist_d, count);
    		MLDouble mlDouble5 = new MLDouble("rssiimpinjlist_d", rdlist_d, count);
    		MLDouble mlDouble6 = new MLDouble("phasedeglist", plist_d, count);
    		MLDouble mlDouble7 = new MLDouble("msgfreqlist", flist_d, count);
    		ArrayList list = new ArrayList();
    		list.add(mlDouble1);
    		list.add(mlDouble2);
    		list.add(mlDouble3);
    		list.add(mlDouble4);
    		list.add(mlDouble5);
    		list.add(mlDouble6);
    		list.add(mlDouble7);
    		String data_name = new String("C:\\Users\\Vicon-OEM\\Desktop\\Kan_Group_folder\\IRB_20191205\\CollectedData\\" + usecase + "_Case_" + caseindex + ".mat");
    		new MatFileWriter(data_name, list);
    		System.out.println("Write to .mat file success.");

    		System.out.println("s1 is: " + s1);
    		System.out.println("s2 is: " + s2);

    		System.out.println("Sleep for 1 second ...");
    		try {
    			Thread.sleep(1000);
    		} catch (InterruptedException ex) {
    			logger.error("Sleep Interruped");
    		} 

        	cnt ++;
    	}

    	masterpiece.disconnect();
    	System.out.println("epclist size is: " + epclist.size()); 



    	// write to excel file
    	// create blank wossrkbook
    	XSSFWorkbook workbook = new XSSFWorkbook();
    	XSSFRow row;
    	Cell cell;

    	// create a blank sheet
    	XSSFSheet spreadsheet = workbook.createSheet("Sheet1");

    	// first write the title of the columns
    	row = spreadsheet.createRow(0);
    	cell = row.createCell(0);
    	cell.setCellValue("Channel Index");
    	cell = row.createCell(1);
    	cell.setCellValue("Tag Index");
    	cell = row.createCell(2);
    	cell.setCellValue("Antenna ID");
    	cell = row.createCell(3);
    	cell.setCellValue("Tag RSSI/dBm");
    	cell = row.createCell(4);
    	cell.setCellValue("Tag RSSI/nW");
    	cell = row.createCell(5);
    	cell.setCellValue("Tag Phase");

    	// then write the received tag data
    	int rowid = 1;
    	while (rowid < epclist.size() + 1) {
    		row = spreadsheet.createRow(rowid);
    		Double rowfreq = ht.get(chindlist.get(rowid - 1).intValue() - 1);
    		int rowfreq_int = rowfreq.intValue();
    		cell = row.createCell(0);
    		cell.setCellValue(chindlist.get(rowid - 1).toString());
    		cell = row.createCell(1);
    		cell.setCellValue(tagindexlist.get(rowid - 1).toString());
    		cell = row.createCell(2);
    		cell.setCellValue(antennalist.get(rowid - 1).toString());
    		cell = row.createCell(3);
    		cell.setCellValue(rssiimpinjlist.get(rowid - 1).toString());
    		cell = row.createCell(4);
    		cell.setCellValue(rssiimpinjlist_d.get(rowid - 1).toString());
    		cell = row.createCell(5);
    		cell.setCellValue(phasedeglist.get(rowid - 1).toString());
    		cell = row.createCell(6);
    		cell.setCellValue(epclist.get(rowid - 1));
    		cell = row.createCell(7);
    		cell.setCellValue(msgfreqlist.get(rowid - 1));
    		rowid += 1;
    	}

    	XSSFSheet spreadsheet2 = workbook.createSheet(" Next Channel List ");
    	row = spreadsheet2.createRow(0);
    	cell = row.createCell(0);
    	cell.setCellValue("Next Channel Index");
    	rowid = 1;
    	while (rowid < nextlist.size() + 1) {
    		row = spreadsheet2.createRow(rowid);
    		cell = row.createCell(0);
    		cell.setCellValue(nextlist.get(rowid - 1).toString());
    		rowid += 1;
    	}

    	XSSFSheet spreadsheet3 = workbook.createSheet(" Hop Table ");
    	row = spreadsheet3.createRow(0);
    	cell = row.createCell(0);
    	cell.setCellValue("Hop Table Index");
    	cell = row.createCell(1);
    	cell.setCellValue("Hop Table Frequency");
    	rowid = 1;
    	while (rowid < ht.size() + 1) {
    		row = spreadsheet3.createRow(rowid);
    		cell = row.createCell(0);
    		cell.setCellValue(hti.get(rowid - 1).toString());
    		cell = row.createCell(1);
    		cell.setCellValue(ht.get(rowid - 1).toString());
    		rowid += 1;
    	}

    	// Write the workbook in file system (with this, we don't need to write to Excel files as specified above)
    	// Codes below output exactly the same data as codes above which write data to Excel.
    	FileOutputStream out = new FileOutputStream(new File("C:\\Users\\Vicon-OEM\\Desktop\\Kan_Group_folder\\IRB_20191205\\CollectedData\\" + usecase + "_Case_" + caseindex + ".xlsx"));
    	workbook.write(out);
    	System.out.println("result_data.xlsx written successfully"); 
    	
    	System.exit(0);
    }
}

