using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Runtime.InteropServices;
using System.Windows.Forms;
namespace MatLib
{
    public class ScreenCapture
    {
        bool IsbuttonCapture;
        //
        private const int WH_KEYBOARD_LL = 13;
        //private const int WH_KEYBOARD_LL = 13;  
        private const int VK_F1 = 0x70;
        //private static LowLevelKeyboardProc _proc = HookCallback;
        private static IntPtr _hookID = IntPtr.Zero;
        public ScreenCapture(bool IsbuttonCapture)
        {
            this.IsbuttonCapture = IsbuttonCapture;
        }
    }
}
