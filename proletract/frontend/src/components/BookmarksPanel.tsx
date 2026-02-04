import React, { useState, useEffect } from 'react';
import { useSession } from '../contexts/SessionContext';
import './BookmarksPanel.css';

export interface Bookmark {
  id: string;
  region: string;
  name: string;
  folder?: string;
  timestamp: number;
  notes?: string;
}

interface BookmarksPanelProps {
  currentRegion?: string;
  onBookmarkSelect?: (region: string) => void;
  onOpenSessionManager?: () => void;
  getSessionData?: () => { vcfPath: string; selectedRegion: string; selectedGenotypes: string[]; publicVcfFolder: string; cohortFolder: string; mode: string };
}

const STORAGE_KEY = 'proletract_bookmarks';

const BookmarksPanel: React.FC<BookmarksPanelProps> = ({ currentRegion, onBookmarkSelect, onOpenSessionManager, getSessionData }) => {
  const [bookmarks, setBookmarks] = useState<Bookmark[]>([]);
  const [isOpen, setIsOpen] = useState(false);
  const [folders, setFolders] = useState<string[]>([]);
  const [selectedFolder, setSelectedFolder] = useState<string>('');
  
  // Session management
  const {
    currentSession,
    saveSession,
    sessions,
  } = useSession();
  
  const [showSessionTooltip, setShowSessionTooltip] = useState(false);
  const sessionTooltipTimerRef = React.useRef<ReturnType<typeof setTimeout> | null>(null);

  useEffect(() => {
    try {
      const stored = localStorage.getItem(STORAGE_KEY);
      if (stored) {
        const allBookmarks: Bookmark[] = JSON.parse(stored);
        setBookmarks(allBookmarks);
        
        const uniqueFolders = Array.from(new Set(allBookmarks.map(b => b.folder).filter(Boolean))) as string[];
        setFolders(uniqueFolders);
      }
    } catch (error) {
      console.error('Error loading bookmarks:', error);
    }
  }, []);

  useEffect(() => {
    return () => {
      if (sessionTooltipTimerRef.current) {
        clearTimeout(sessionTooltipTimerRef.current);
      }
    };
  }, []);

  const saveBookmarks = (updatedBookmarks: Bookmark[]) => {
    try {
      localStorage.setItem(STORAGE_KEY, JSON.stringify(updatedBookmarks));
      setBookmarks(updatedBookmarks);
    } catch (error) {
      console.error('Error saving bookmarks:', error);
    }
  };

  const handleAddBookmark = () => {
    if (!currentRegion) {
      alert('No region selected to bookmark');
      return;
    }

    const name = prompt('Enter bookmark name:', currentRegion);
    if (!name) return;

    const folder = selectedFolder || undefined;
    const newBookmark: Bookmark = {
      id: `bookmark_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
      region: currentRegion,
      name: name.trim(),
      folder,
      timestamp: Date.now(),
    };

    const updated = [...bookmarks, newBookmark];
    saveBookmarks(updated);
    
    if (folder && !folders.includes(folder)) {
      setFolders(prev => [...prev, folder]);
    }
  };

  const handleDeleteBookmark = (id: string) => {
    if (window.confirm('Delete this bookmark?')) {
      const updated = bookmarks.filter(b => b.id !== id);
      saveBookmarks(updated);
    }
  };

  const handleSelectBookmark = (bookmark: Bookmark) => {
    if (onBookmarkSelect) {
      onBookmarkSelect(bookmark.region);
    }
    setIsOpen(false);
  };

  const filteredBookmarks = selectedFolder
    ? bookmarks.filter(b => b.folder === selectedFolder)
    : bookmarks;

  const bookmarksByFolder = folders.reduce((acc, folder) => {
    acc[folder] = bookmarks.filter(b => b.folder === folder);
    return acc;
  }, {} as Record<string, Bookmark[]>);

  const unorganizedBookmarks = bookmarks.filter(b => !b.folder);

  const handleSaveSession = (e: React.MouseEvent) => {
    // If Ctrl/Cmd is held, open full session manager
    if (e.ctrlKey || e.metaKey) {
      if (onOpenSessionManager) {
        onOpenSessionManager();
      }
      return;
    }
    
    if (!currentSession) {
      const sessionName = prompt('Enter session name:');
      if (!sessionName) return;
      
      const data = getSessionData?.() ?? {
        vcfPath: '',
        selectedRegion: currentRegion || '',
        selectedGenotypes: [] as string[],
        publicVcfFolder: '',
        cohortFolder: '',
        mode: 'individual' as const,
      };
      saveSession({
        name: sessionName,
        vcfPath: data.vcfPath,
        selectedRegion: data.selectedRegion,
        selectedGenotypes: data.selectedGenotypes,
        publicVcfFolder: data.publicVcfFolder,
        cohortFolder: data.cohortFolder,
        mode: data.mode === 'cohort-read' || data.mode === 'cohort-assembly' ? data.mode : 'individual',
      });
    } else {
      saveSession({
        name: currentSession.name,
        vcfPath: currentSession.vcfPath,
        selectedRegion: currentSession.selectedRegion,
        selectedGenotypes: currentSession.selectedGenotypes,
        publicVcfFolder: currentSession.publicVcfFolder,
        cohortFolder: currentSession.cohortFolder,
        mode: currentSession.mode,
        notes: currentSession.notes,
      });
    }
  };

  const formatDate = (timestamp: number) => {
    return new Date(timestamp).toLocaleString();
  };

  const getSessionTooltipText = () => {
    if (!currentSession) {
      return 'No active session\nClick to save current state\nCtrl+Click to open session manager';
    }
    return `Session: ${currentSession.name}\nSaved: ${formatDate(currentSession.timestamp)}\n${sessions.length} total sessions\nCtrl+Click to open session manager`;
  };

  return (
    <div className="bookmarks-panel">
      <div className="bookmarks-actions-group">
        <button
          className="bookmarks-toggle-button"
          onClick={() => setIsOpen(!isOpen)}
          title="Bookmarks"
        >
          <span className="bookmarks-icon">🔖</span>
          <span className="bookmarks-label">Bookmarks</span>
          {bookmarks.length > 0 && (
            <span className="bookmarks-count">{bookmarks.length}</span>
          )}
        </button>
        
        <div 
          className="session-save-button-wrapper"
          onMouseEnter={() => {
            sessionTooltipTimerRef.current = setTimeout(() => setShowSessionTooltip(true), 400);
          }}
          onMouseLeave={() => {
            if (sessionTooltipTimerRef.current) {
              clearTimeout(sessionTooltipTimerRef.current);
              sessionTooltipTimerRef.current = null;
            }
            setShowSessionTooltip(false);
          }}
        >
          <button
            type="button"
            className="session-save-button"
            onClick={handleSaveSession}
            title={getSessionTooltipText()}
          >
            <span className="session-save-icon">💾</span>
          </button>
          {showSessionTooltip && (
            <div className="session-save-tooltip">
              <div className="session-tooltip-content">
                {currentSession ? (
                  <>
                    <div className="session-tooltip-line"><strong>Session:</strong> {currentSession.name}</div>
                    <div className="session-tooltip-line"><strong>Saved:</strong> {formatDate(currentSession.timestamp)}</div>
                    {currentSession.selectedRegion && (
                      <div className="session-tooltip-line"><strong>Region:</strong> {currentSession.selectedRegion}</div>
                    )}
                    <div className="session-tooltip-line"><strong>Total Sessions:</strong> {sessions.length}</div>
                  </>
                ) : (
                  <div className="session-tooltip-line">No active session</div>
                )}
                <div className="session-tooltip-line session-tooltip-hint">Click to save current state</div>
              </div>
            </div>
          )}
        </div>
      </div>

      {isOpen && (
        <div className="bookmarks-dialog">
          <div className="bookmarks-dialog-header">
            <h3>🔖 Bookmarks</h3>
            <button
              className="bookmarks-close-button"
              onClick={() => setIsOpen(false)}
            >
              ✕
            </button>
          </div>

          <div className="bookmarks-dialog-content">
            {currentRegion && (
              <button
                className="bookmarks-add-button"
                onClick={handleAddBookmark}
              >
                + Bookmark Current Region
              </button>
            )}

            {folders.length > 0 && (
              <div className="bookmarks-folders">
                <label>Filter by folder:</label>
                <select
                  className="bookmarks-folder-select"
                  value={selectedFolder}
                  onChange={(e) => setSelectedFolder(e.target.value)}
                >
                  <option value="">All Bookmarks</option>
                  {folders.map(folder => (
                    <option key={folder} value={folder}>{folder}</option>
                  ))}
                </select>
              </div>
            )}

            {filteredBookmarks.length === 0 ? (
              <div className="bookmarks-empty">
                {bookmarks.length === 0
                  ? 'No bookmarks yet'
                  : 'No bookmarks in this folder'}
              </div>
            ) : (
              <div className="bookmarks-list">
                {selectedFolder ? (
                  // Show bookmarks in selected folder
                  bookmarksByFolder[selectedFolder]?.map(bookmark => (
                    <div key={bookmark.id} className="bookmark-item">
                      <div
                        className="bookmark-item-content"
                        onClick={() => handleSelectBookmark(bookmark)}
                      >
                        <div className="bookmark-item-name">{bookmark.name}</div>
                        <div className="bookmark-item-region">{bookmark.region}</div>
                        <div className="bookmark-item-time">
                          {new Date(bookmark.timestamp).toLocaleString()}
                        </div>
                      </div>
                      <button
                        className="bookmark-item-delete"
                        onClick={() => handleDeleteBookmark(bookmark.id)}
                        title="Delete bookmark"
                      >
                        🗑️
                      </button>
                    </div>
                  ))
                ) : (
                  // Show all bookmarks organized by folder
                  <>
                    {folders.map(folder => (
                      <div key={folder} className="bookmark-folder-section">
                        <div className="bookmark-folder-header">{folder}</div>
                        {bookmarksByFolder[folder]?.map(bookmark => (
                          <div key={bookmark.id} className="bookmark-item">
                            <div
                              className="bookmark-item-content"
                              onClick={() => handleSelectBookmark(bookmark)}
                            >
                              <div className="bookmark-item-name">{bookmark.name}</div>
                              <div className="bookmark-item-region">{bookmark.region}</div>
                            </div>
                            <button
                              className="bookmark-item-delete"
                              onClick={() => handleDeleteBookmark(bookmark.id)}
                            >
                              🗑️
                            </button>
                          </div>
                        ))}
                      </div>
                    ))}
                    {unorganizedBookmarks.length > 0 && (
                      <div className="bookmark-folder-section">
                        <div className="bookmark-folder-header">Unorganized</div>
                        {unorganizedBookmarks.map(bookmark => (
                          <div key={bookmark.id} className="bookmark-item">
                            <div
                              className="bookmark-item-content"
                              onClick={() => handleSelectBookmark(bookmark)}
                            >
                              <div className="bookmark-item-name">{bookmark.name}</div>
                              <div className="bookmark-item-region">{bookmark.region}</div>
                            </div>
                            <button
                              className="bookmark-item-delete"
                              onClick={() => handleDeleteBookmark(bookmark.id)}
                            >
                              🗑️
                            </button>
                          </div>
                        ))}
                      </div>
                    )}
                  </>
                )}
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  );
};

export default BookmarksPanel;




