from Tkinter import Tk, Label, Button

class DebateGUI:

    def __init__(self, master):

        self.master = master

        master.title("Debate Flow")

        self.label = Label(master, text="Choose One:")
        self.label.pack()

        self.con_button = Button(master, text="Contention", command=self.contention)
        self.con_button.pack()

        self.close_button = Button(master, text="Close", command=master.quit)
        self.close_button.pack()

        
    def contention(self):
        self.newcon_button = Button(self, text="What's the contention", command=self.newcon)
    def newcon(self):
        print(stupid)




root = Tk()
my_gui = DebateGUI(root)
root.mainloop()