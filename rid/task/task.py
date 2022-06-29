import os
from collections.abc import (
    Sequence,
)
from typing import (
    List,
    Tuple,
    Dict,
    Optional,
    Any
)

class Task():
    """Define the files needed by an MD task. 
    Examples
    --------
    >>> # this example dumps all files needed by the task.
    >>> files = exploration_task.files()
    ... for file_name, file_content in files.items():
    ...     with open(file_name, 'w') as fp:
    ...         fp.write(file_content)    
    """
    def __init__(
            self, 
    ):
        self._files = {}
    
    def add_property(
            self,
            properties: Dict
        ):
        self.__dict__.update(properties)

    def add_file(
            self,
            fname : str,
            fcont : Tuple[Optional[str, bytes], str],
    ):
        """Add file to the task
        
        Parameters
        ----------
        fname : str
            The name of the file
        fcont : str
            The content of the file.
        """
        self._files[fname] = fcont
        return self
    
    @property
    def files(self) -> Dict:
        """Get all files for the task.
        
        Returns
        -------
        files : dict
            The dict storing all files for the task. The file name is a key of the dict, and the file content is the corresponding value.
        """
        return self._files


class TaskGroup(Sequence):
    """A group of exploration tasks. Implemented as a `list` of `MDTask`.
    """
    def __init__(self):
        super().__init__()
        self.clear()

    def __getitem__(self, ii:int) -> Task:
        """Get the `ii`th task"""
        return self.task_list[ii]

    def __len__(self) -> int:
        """Get the number of tasks in the group"""
        return len(self.task_list)

    def clear(self)->None:
        self._task_list = []

    @property
    def task_list(self) -> List[Task]:
        """Get the `list` of `MDTask`""" 
        return self._task_list

    def add_task(self, task: MDTask):
        """Add one task to the group."""
        self.task_list.append(task)
        return self

    def add_group(
            self,
            group : 'MDTaskGroup',
    ):
        """Add another group to the group."""
        self._task_list = self._task_list + group._task_list
        return self

    def __add__(
            self,
            group : 'MDTaskGroup',
    ):        
        """Add another group to the group."""
        return self.add_group(group)
