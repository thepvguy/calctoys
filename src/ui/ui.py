import tkinter as tk


def tk_sticky_all():
    return tk.N+tk.E+tk.W+tk.S


class MainWindow(tk.Frame):
    def __init__(self, parent):
        super(MainWindow, self).__init__()
        self.parent = parent
        self.component_choice_pane = None

        self.__MainWindowSetup()

    def __MainWindowSetup(self):
        self.component_choice_pane = ComponentChoicePane(self)
        self.component_choice_pane.grid(row=0, column=0, sticky=tk_sticky_all())


class ComponentChoicePane(tk.LabelFrame):
    def __init__(self, parent):
        tk.LabelFrame.__init__(self, parent)
        self.parent = parent
        self.text = "Component Choices"
        self.thing = tk.Label(text="wubba lubba dub dub").pack()


def init_ui(context=None):
    root = tk.Tk()
    MainWindow(root).grid(row=0, column=0, sticky=tk_sticky_all())
    root.mainloop()
