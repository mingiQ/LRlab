try:
    import pyvisa

except Exception as e:
    print(e)
    print("Warning VISA library import failed")
import telnetlib
import socket
import time

try:
    import serial
except ImportError:
    print("Warning serial library import failed.")


class Instrument(object):
    """
    A subclass of Instrument is an instrument which communicates over a certain
    channel. The subclass must define the methods write and read, for
    communication over that channel
    """
    address = ''  # Address of instrument
    name = ''  # Instrument Name
    enabled = False  # If enabled=False commands should not be sent
    instrument_type = ''  # Instrument type
    protocol = ''  # Protocol
    id_string = ''  # id string
    query_sleep = 0  # seconds to wait between write and read
    term_char = '\n'  # character to be appended to all writes

    # operation_range={}        #map to hold the operation range

    def __init__(self, name, address='', enabled=True, timeout=1, query_sleep=0):
        """
        :param name:
        :param address:
        :param enabled:
        :param timeout: timeout for low-level queries in seconds
        :return:
        """
        self.name = name
        self.address = address
        self.enabled = enabled
        self.timeout = timeout  # timeout for connection, different from timeout for query
        self.query_sleep = query_sleep

    def get_name(self):
        return self.name

    def get_id(self):
        return "Default Instrument %s" % (self.name)

    def encode_s(self, s):
        if type(self.term_char) == str:
            term_char = str.encode(self.term_char)
        else:
            term_char = self.term_char

        if type(s) == str:
            return str.encode(s) + term_char
        else:
            return s + term_char

    def query(self, cmd, timeout=None):
        self.write(cmd)
        time.sleep(self.query_sleep)
        return self.read(timeout)

    def queryb(self, cmd, timeout=None):
        self.write(cmd)
        time.sleep(self.query_sleep)
        return self.readb(timeout)

    def set_timeout(self, timeout=None):
        if timeout is not None:
            self.timeout = timeout

    def get_timeout(self):
        return self.timeout

    def set_query_sleep(self, query_sleep):
        self.query_sleep = query_sleep

    def get_query_sleep(self):
        return self.query_sleep

    def get_settings(self):
        settings = {}
        settings['name'] = self.name
        settings['address'] = self.address
        settings['instrument_type'] = self.instrument_type
        settings['protocol'] = self.protocol
        return settings

    def set_settings(self, settings):
        print(settings)

    def attr(self, name):
        "re-naming of __getattr__ which is unavailable when proxied"
        return getattr(self, name)

# =============================================================================
# # In[]
# inst = Instrument(name='test')
# 
# print(inst.encode_s('*IDN?'))
# 
# # In[]
# =============================================================================
class VisaInstrument(Instrument):
    def __init__(self, name, address='', enabled=True, timeout=1.0, rm=False, **kwargs):
        Instrument.__init__(self, name, address, enabled, timeout, **kwargs)
        if self.enabled:
            self.protocol = 'VISA'
            self.timeout = timeout
            address = address.upper()
            if rm is not False:
                self.instrument = pyvisa.ResourceManager(rm).open_resource(address)
            else:
                self.instrument = pyvisa.ResourceManager().open_resource(address)
            self.instrument.timeout = timeout * 1000

    def _query(self, s):
        if self.enabled: return(self.instrument.query(self.encode_s(s)))    

    def write(self, s):
        if self.enabled: self.instrument.write(self.encode_s(s))

    def read(self, timeout=None):
        # todo: implement timeout, reference SocketInstrument.read
        if self.enabled: return self.instrument.read().decode()

    def readb(self, timeout=None):
        # todo: implement timeout, reference SocketInstrument.read
        if self.enabled: return self.instrument.read()

    def close(self):
        if self.enabled: self.instrument.close()
        
    def read_lineb(self, eof_char=b'\n', timeout=None):
        done = False
        while done is False:
            buffer_str = self.readb(timeout)
            # print "buffer_str", [buffer_str]
            if buffer_str is None:
                pass  # done = True
            elif buffer_str[-len(eof_char):] == eof_char:
                done = True
                yield buffer_str
            else:
                yield buffer_str


class TelnetInstrument(Instrument):
    def __init__(self, name, address='', enabled=True, timeout=10, **kwargs):
        Instrument.__init__(self, name, address, enabled, timeout, **kwargs)
        self.protocol = 'Telnet'
        if len(address.split(':')) > 1:
            self.port = int(address.split(':')[1])
        if self.enabled:
            self.tn = telnetlib.Telnet(address.split(':')[0], self.port)

    def write(self, s):
        if self.enabled: self.tn.write(self.encode_s(s))

    def read(self, timeout=None):
        # todo: implement timeout, reference SocketInstrument.read
        if self.enabled: return self.tn.read_some().decode()

    def readb(self, timeout=None):
        # todo: implement timeout, reference SocketInstrument.read
        if self.enabled: return self.tn.read_some()

    def close(self):
        if self.enabled: self.tn.close()


import select


class SocketInstrument(Instrument):
    default_port = 23

    def __init__(self, name, address='', enabled=True, recv_length=1024, timeout=1.0, **kwargs):
        Instrument.__init__(self, name, address, enabled, timeout, **kwargs)
        self.protocol = 'socket'
        self.recv_length = recv_length
        if len(address.split(':')) > 1:
            self.port = int(address.split(':')[1])
            self.ip = address.split(':')[0]
        else:
            self.ip = address
            self.port = self.default_port
        self.on_enable()

    def on_enable(self):
        if self.enabled:
            self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.socket.connect((self.ip, self.port))
            self.set_timeout(self.timeout)
            self.socket.setblocking(0)

    def set_enable(self, enable=True):
        self.enabled = enable
        self.on_enable()

    def set_timeout(self, timeout):
        Instrument.set_timeout(self, timeout)
        if self.enabled: self.socket.settimeout(self.timeout)

    def write(self, s):
        self.cmd = str(s)
        if self.enabled: self.socket.send(self.cmd.encode())

    # def query(self, s):
    #     self.write(s)
    #     time.sleep(self.query_sleep)
    #     return self.read()

    def read(self, timeout=None):
        if timeout == None: timeout = self.timeout
        ready = select.select([self.socket], [], [], timeout)
        if (ready[0] and self.enabled):
            return self.socket.recv(self.recv_length).decode()

    def readb(self, timeout=None):
        if timeout == None: timeout = self.timeout
        ready = select.select([self.socket], [], [], timeout)
        if (ready[0] and self.enabled):
            return self.socket.recv(self.recv_length)

    def read_line(self, eof_char=b'\n', timeout=None):
        done = False
        while done is False:
            buffer_str = self.read(timeout)
            # print "buffer_str", [buffer_str]
            if buffer_str is None:
                pass  # done = True
            elif buffer_str[-len(eof_char):] == eof_char:
                done = True
                yield buffer_str
            else:
                yield buffer_str

    def read_lineb(self, eof_char=b'\n', timeout=None):
        done = False
        while done is False:
            buffer_str = self.readb(timeout)
            # print "buffer_str", [buffer_str]
            if buffer_str is None:
                pass  # done = True
            elif buffer_str[-len(eof_char):] == eof_char:
                done = True
                yield buffer_str
            else:
                yield buffer_str


class SerialInstrument(Instrument):
    # todo: the `baudrate` and `querysleep` need to be updated to band_rate and query_sleep
    def __init__(self, name, address, enabled=True, timeout=1.0,
                 recv_length=1024, baudrate=9600, query_sleep=1.0):
        Instrument.__init__(self, name, address, enabled)
        self.protocol = 'serial'
        self.enabled = enabled
        if self.enabled:
            try:
                self.ser = serial.Serial(port=address, baudrate=baudrate, timeout=timeout)
            except serial.SerialException:
                print('Cannot create a connection to port ' + str(address) + '.\n')
        self.set_timeout(timeout)
        self.recv_length = recv_length
        self.query_sleep = query_sleep

    def set_timeout(self, timeout):
        Instrument.set_timeout(self, timeout)
        if self.enabled: self.ser.timeout = self.timeout

    def test(self):
        self.ser.setTimeout(self.timeout)

    def write(self, s):
        cmd = self.encode_s(s)
        if self.enabled: self.ser.write(cmd)
        
    def writeb(self, s):
        cmd = s+'\n'
        self.ser.write(str.encode(cmd))

    def read(self, timeout=None):
        # todo: implement timeout, reference SocketInstrument.read
        if self.enabled: return self.ser.read(self.recv_length).decode()

    def readb(self, timeout=None):
        # todo: implement timeout, reference SocketInstrument.read
        if self.enabled: 
            data = self.ser.readline().rstrip()
        return data.decode()

    def reset_connection(self):
        self.ser.close()
        time.sleep(self.query_sleep)
        self.ser.open()

    def flush(self):
        self.ser.flush()
        time.sleep(1)

    def flush_input(self):
        self.ser.flushInput()
        time.sleep(1)

    def flush_output(self):
        self.ser.flushOutput()
        time.sleep(1)

    def __del__(self):
        try:
            self.ser.close()
        except Exception as e:
            print(e)
            print('cannot properly close the serial connection.')